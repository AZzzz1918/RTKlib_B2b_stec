/*------------------------------------------------------------------------------
* solution.c : solution functions
*-----------------------------------------------------------------------------*/
#include <ctype.h>
#include "gamp.h"

/* constants and macros ------------------------------------------------------*/

/* solution option to field separator ----------------------------------------*/
static const char *opt2sep(const solopt_t *opt)
{
	if (!*opt->sep) return " ";
	else if (!strcmp(opt->sep,"\\t")) return "\t";
	return opt->sep;
}
/* output solution as the form of lat/lon/height, modified by fzhou @ GFZ, 2017-01-23 -----------------------------*/
static int outpos(unsigned char *buff, const char *s, const sol_t *sol, const solopt_t *opt)
{
	double pos[3],dxyz[3],denu[3];
	const char *sep=opt2sep(opt);
	char *p=(char *)buff;
	int i;
	double ep[6];
	gtime_t time;
	int week;
	double sow;

	ep[0]=str2num(s,0,4);
	ep[1]=str2num(s,5,2);
	ep[2]=str2num(s,8,2);
	ep[3]=str2num(s,11,2);
	ep[4]=str2num(s,14,2);
	ep[5]=str2num(s,17,2);
	time=epoch2time(ep);
	sow=time2gpst(time,&week);

	for (i=0;i<3;i++) dxyz[i]=denu[i]=0.0;

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s%14.4f%s%14.4f%s%14.4f",
		(int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep,
		sol->rr[0],sep,sol->rr[1],sep,sol->rr[2]);

	if (PPP_Glo.crdTrue[0]==0.0) ecef2pos(sol->rr,pos);
	else ecef2pos(PPP_Glo.crdTrue,pos);

	if (PPP_Glo.crdTrue[0]==0.0) 
		denu[0]=denu[1]=denu[2]=0.0;
	else {
		ecef2pos(PPP_Glo.crdTrue,pos);
		for (i=0;i<3;i++) 
			denu[i]=dxyz[i]=sol->rr[i]-PPP_Glo.crdTrue[i];
		ecef2enu(pos,dxyz,denu);
	}

	p+=sprintf(p,"%s%8.4f%s%8.4f%s%8.4f%s%8.4f",sep,denu[0],sep,denu[1],sep,denu[2],
		sep,norm(denu,3));

	p+=sprintf(p,"\n");

	return p-(char *)buff;
}
/* output solution body --------------------------------------------------------
* output solution body to buffer
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
static int outsols(unsigned char *buff, const sol_t *sol, const solopt_t *opt)
{
    gtime_t time;
    double gpst;
    int week,timeu;
    const char *sep=opt2sep(opt);
    char s[64];
    unsigned char *p=buff;

    if (PPP_Glo.prcOpt_Ex.posMode>=PMODE_PPP_KINEMA&&sol->stat==SOLQ_SINGLE) {
        p+=sprintf((char *)p,"\n");
        return p-buff;
    }
    
    if (sol->stat<=SOLQ_NONE) {
        p+=sprintf((char *)p,"\n");
        return p-buff;
    }
    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);
    
    time=sol->time;
    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    
    if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }

	p+=outpos(p,s,sol,opt);

    return p-buff;
}
/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsol(FILE *fp, const sol_t *sol, const solopt_t *opt, int isol)
{
    unsigned char buff[MAXSOLMSG+1]={'\0'};
    int n;
    
	if ((n=outsols(buff,sol,opt))>0) {
		fwrite(buff,n,1,fp);
    }
}

//output Debug information
extern void outDebug(int bWinOut, int bFpOut, int bTMOut)
{
	int bOk=0;

	if (bWinOut) {
		if (bTMOut) printf("%s  ",PPP_Glo.chTime);
		printf("%s",PPP_Glo.chMsg);
	}

	if (PPP_Glo.outFp[OFILE_DEBUG]==NULL) return;

	if (PPP_Glo.prcOpt_Ex.solType<=1)
		bOk=1;
	else if (PPP_Glo.prcOpt_Ex.solType==2) {
		if (PPP_Glo.revs==0) bOk=1;
	}
	else if (PPP_Glo.prcOpt_Ex.solType==3) {
		if (PPP_Glo.revs==1) bOk=1;
	}
	else if (PPP_Glo.prcOpt_Ex.solType==4) {
		if ((PPP_Glo.revs==0)&&(PPP_Glo.iTag==2)) bOk=1;
	}

	if (bOk&&bFpOut) {
		if (bTMOut) fprintf(PPP_Glo.outFp[OFILE_DEBUG],"%s  ",PPP_Glo.chTime);
		fprintf(PPP_Glo.outFp[OFILE_DEBUG],"%s",PPP_Glo.chMsg);
	}
}

//output residuals
static void outResi(FILE *fp, rtk_t *rtk, gtime_t time, int frq, int itype, int yanqian)
{
    unsigned char buff[MAXSOLMSG+1];
    int i,j,n,i0,i1;
    char *p=(char *)buff;
    char *sep = " ";
	double sow,ep[6];
	int week;
	double resc_pri[NFREQ],resp_pri[NFREQ],resc_pos[NFREQ],resp_pos[NFREQ];

	for (i=0;i<3;i++) {
		resc_pri[i]=0.0;
		resp_pri[i]=0.0;
		resc_pos[i]=0.0;
		resp_pos[i]=0.0;
	}

	time2epoch(time,ep);
	sow=time2gpst(time,&week);
    
	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],sep,
		(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);
    
	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		if (rtk->ssat[i].vsat[0]==1) {
			for (j=0;j<NFREQ;j++) {
				resc_pri[j]=rtk->ssat[i].resc_pri[j];
				resp_pri[j]=rtk->ssat[i].resp_pri[j];
				resc_pos[j]=rtk->ssat[i].resc_pos[j];
				resp_pos[j]=rtk->ssat[i].resp_pos[j];
			}
		}
		else if (rtk->ssat[i].vsat[0]==0) {
			for (j=0;j<NFREQ;j++) {
				resc_pri[j]=99999.0;
				resp_pri[j]=99999.0;
				resc_pos[j]=99999.0;
				resp_pos[j]=99999.0;
			}
		}
		if (yanqian) {
			if (frq==1) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pri[0],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pri[0],sep);
			}
			else if (frq==2) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pri[1],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pri[1],sep);
			}
			else if (frq==3) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pri[2],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pri[2],sep);
			}
		}
		else {
			if (frq==1) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pos[0],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pos[0],sep);
			}
			else if (frq==2) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pos[1],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pos[1],sep);
			}
			else if (frq==3) {
				if (itype==0) p+=sprintf(p,"%10.4f%s",resc_pos[2],sep);
				if (itype==1) p+=sprintf(p,"%10.4f%s",resp_pos[2],sep);
			}
		}
	}

    p+=sprintf(p,"\n");
 
    n=p-(char *)buff;
    
    if (n>0) {
        fwrite(buff,n,1,fp);
    }
}

//output LC minus PC observations
static void outLCmPC(FILE *fp, rtk_t *rtk, gtime_t time)
{
    unsigned char buff[MAXSOLMSG+1];
    int i,n,i0,i1;
    char *p=(char *)buff;
    char *sep = " ";
	double sow,ep[6];
	int week;
	double lcmpc=0.0;

	time2epoch(time,ep);
	sow=time2gpst(time,&week);
    
	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],sep,
		(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);
    
	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		if (rtk->ssat[i].vsat[0]==1) {
            lcmpc = rtk->ssat[i].LC-rtk->ssat[i].PC;
		}
		else if (rtk->ssat[i].vsat[0]==0) {
            lcmpc = 99999.0;
		}
        p+=sprintf(p,"%10.4f%s",lcmpc,sep);
	}

    p+=sprintf(p,"\n");
 
    n=p-(char *)buff;
    
    if (n>0) {
        fwrite(buff,n,1,fp);
    }
}

//output number of satellite and PDOP
static void outPdop(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int n;
	char *p=(char *)buff;
	char *sep = " ";
	double sow,ep[6];
	int week,nsat;

	time2epoch(time,ep);

	sow=time2gpst(time,&week);

	nsat=rtk->sol.ns[0]+rtk->sol.ns[1]+rtk->sol.ns[2]+rtk->sol.ns[3]+rtk->sol.ns[4];
	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s%3d%s%2d%s%2d%s%2d%s%2d%s%2d%s%10.3f",
		(int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,
		sep,sow,sep,nsat,sep,rtk->sol.ns[0],sep,rtk->sol.ns[1],sep,rtk->sol.ns[2],sep,rtk->sol.ns[3],sep,rtk->sol.ns[4],
		sep,rtk->sol.dop[1]);

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

//output inter-system biases (ISBs) for multi-GNSS
static void outIsb(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int j=0,n;
	char *p=(char *)buff;
	char *sep = " ";
	double sow,ep[6],isb;
	int week;

	time2epoch(time,ep);

	sow=time2gpst(time,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],
		sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	int gps_en = (rtk->opt.navsys & 0x1) ? 1 : 0;
	int glo_en = (rtk->opt.navsys >> 2 & 0x1) ? 1 : 0;
	int gal_en = (rtk->opt.navsys >> 3 & 0x1) ? 1 : 0;
	int bds_en = (rtk->opt.navsys >> 5 & 0x1) ? 1 : 0;
	int qzs_en = (rtk->opt.navsys >> 4 & 0x1) ? 1 : 0;

	if (gps_en && glo_en) {
		j = IC(1, &rtk->opt);
		isb = rtk->x[j] / CLIGHT * 1e+9;
		p += sprintf(p, "%9.3f%s", isb, sep);
	}
	if (gps_en && bds_en) {
		j = IC(2, &rtk->opt);
		isb = rtk->x[j] / CLIGHT * 1e+9;
		p += sprintf(p, "%9.3f%s", isb, sep);
	}
	if (gps_en && gal_en) {
		j = IC(3, &rtk->opt);
		isb = rtk->x[j] / CLIGHT * 1e+9;
		p += sprintf(p, "%9.3f%s", isb, sep);
	}
	if (gps_en && qzs_en) {
		j = IC(4, &rtk->opt);
		isb = rtk->x[j] / CLIGHT * 1e+9;
		p += sprintf(p, "%9.3f%s", isb, sep);
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

//output inter-system biases (ISBs) every 30 min for multi-GNSS
static void outIsb_m(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int n;
	char *p=(char *)buff;
	char *sep = " ";
	double sow,ep[6],isb;
	int week;
	double tdif;

	tdif=timediff(PPP_Glo.t_30min,time);
	if (tdif!=0.0) return;

	time2epoch(time,ep);

	sow=time2gpst(time,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],
		sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	int gps_en = (rtk->opt.navsys & 0x1) ? 1 : 0;
	int glo_en = (rtk->opt.navsys >> 2 & 0x1) ? 1 : 0;
	int gal_en = (rtk->opt.navsys >> 3 & 0x1) ? 1 : 0;
	int bds_en = (rtk->opt.navsys >> 5 & 0x1) ? 1 : 0;
	int qzs_en = (rtk->opt.navsys >> 4 & 0x1) ? 1 : 0;

	if (gps_en && (glo_en || gal_en || bds_en || qzs_en)) {
		isb = PPP_Glo.isb_30min;
		p += sprintf(p, "%9.3f%s", isb, sep);
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

//output inter-frequency biases (IFBs) for GLONASS
static void outIfb(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,j=0,n;
	char *p=(char *)buff;
	char *sep = " ";
	double sow,ep[6],ifb;
	int week;

	time2epoch(time,ep);

	sow=time2gpst(time,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],
		sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	if ((SYS_GLO&rtk->opt.navsys)) {
		n=NICB(&rtk->opt);
		for (i=0;i<n;i++) {
			j=IICB(i+1,&rtk->opt);
			ifb=rtk->x[j];
			p+=sprintf(p,"%9.3f%s",ifb,sep);
		}
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

/* output ZTD */
static void outDtrp(FILE *fp, rtk_t *rtk, gtime_t time)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    char *p=(char *)buff;
    char *sep = " ";
    double sow,ep[6],ztd;
    int week;

    time2epoch(time,ep);

    sow = time2gpst(time,&week);

    p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",
        (int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

    ztd = PPP_Glo.zhd+rtk->x[IT(&rtk->opt)];
    p+=sprintf(p, "%9.4f%s%9.4f%s%9.4f", PPP_Glo.zhd,sep,rtk->x[IT(&rtk->opt)],sep,ztd);

    p+=sprintf(p,"\n");

    n = p-(char *)buff;

    if (n>0) {
        fwrite(buff,n,1,fp);
    }
}

/* output sTEC -------------------*/
static void outStec(FILE* fp, rtk_t *rtk, gtime_t t)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,j=0,n,i0,i1;
	char *p=(char *)buff;
	char *sep=" ";
	double sow,ep[6];
	int week;
	double ion=0.0,deg,mindeg=15.0;
	double fact;

	time2epoch(t,ep);

	sow=time2gpst(t,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],
		sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

    i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }
    
    for (i=i0;i<i1;i++) {
        j=II(i+1,&rtk->opt);
        deg=rtk->ssat[i].azel[1]*R2D;
        if (deg>=mindeg) {
            if (rtk->ssat[i].vsat[0]==1) {
                if (rtk->x[j] == 0.0) {
                    ion = 99999.0;
                }
                else {
                    int sys_tmp = PPP_Glo.sFlag[i].sys;
                    double freq = 0.0;
                    if (sys_tmp == SYS_GLO && PPP_Glo.lam[i][0] != 0.0) {
                        freq = CLIGHT / PPP_Glo.lam[i][0];
                        fact = 40.30E16 / (freq * freq);
                    }
                    else if (sys_tmp == SYS_GPS || sys_tmp == SYS_QZS) {
                        fact = 40.30E16 / FREQ1 / FREQ1;
                    }
                    else if (sys_tmp == SYS_GAL) {
                        fact = 40.30E16 / FREQ1 / FREQ1;
                    }
                    else if (sys_tmp == SYS_CMP) {
                        fact = 40.30E16 / FREQ1_CMP / FREQ1_CMP;
                    }
                    else {
                        fact = 40.30E16 / FREQ1 / FREQ1;
                    }
                    ion = fabs(rtk->x[j] / fact);
                }
            }
            else if (rtk->ssat[i].vsat[0]==0) {
                ion=99999.0;
            }
        }
        else ion=99999.0;
        
        p+=sprintf(p,"%9.3f%s",ion,sep);
    }

	p+=sprintf(p,"\n");

	n = p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}
/* output vTEC -------------------*/
static void outVtec(FILE* fp, rtk_t* rtk, gtime_t t)
{
	unsigned char buff[MAXSOLMSG + 1];
	int i, j = 0, n, i0, i1;
	char* p = (char*)buff;
	char* sep = " ";
	double sow, ep[6];
	int week;
	double ion = 0.0, deg, mindeg = 15.0;
	double elev = 99999.0, vtec = 99999.0, mf = 0.0;
	double Hopt = 506.7, a = 0.9782, R = 6378.137; 
	double pos[3], posp[3] = { 99999.0 }, azel[2]; 
	double cosaz, rp, ap, sinap, tanap;
	extern double Hion, re;
	double fact;

	time2epoch(t, ep);

	sow = time2gpst(t, &week);

	p += sprintf(p, "%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s", (int)ep[0], sep, (int)ep[1],
		sep, (int)ep[2], sep, (int)ep[3], sep, (int)ep[4], sep, (int)ep[5], sep, week, sep, sow, sep);

	i0 = MAXSAT; i1 = 0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i = i0; i < i1; i++) {
		j = II(i + 1, &rtk->opt);
		deg = rtk->ssat[i].azel[1] * R2D;
		if (deg >= mindeg) {
			if (rtk->ssat[i].vsat[0] == 1) {
				if (rtk->x[j] == 0.0) {
					ion = 99999.0;
				}
				else {
					int sys_tmp = PPP_Glo.sFlag[i].sys;
					double freq = 0.0;
					if (sys_tmp == SYS_GLO && PPP_Glo.lam[i][0] != 0.0) {
						freq = CLIGHT / PPP_Glo.lam[i][0];
						fact = 40.30E16 / (freq * freq);
					}
					else if (sys_tmp == SYS_GPS || sys_tmp == SYS_QZS) {
						fact = 40.30E16 / FREQ1 / FREQ1;
					}
					else if (sys_tmp == SYS_GAL) {
						fact = 40.30E16 / FREQ1 / FREQ1;
					}
					else if (sys_tmp == SYS_CMP) {
						fact = 40.30E16 / FREQ1_CMP / FREQ1_CMP;
					}
					else {
						fact = 40.30E16 / FREQ1 / FREQ1;
					}
					ion = fabs(rtk->x[j] / fact);
				}
				elev = rtk->ssat[i].azel[1] * R2D;
				elev = elev * PI / 180.0;
				elev = PI / 2 - elev;
				mf = 1 / (cos(asin(R * sin(a * elev) / (R + Hopt)))); 
				vtec = ion / mf;

				ecef2pos(PPP_Glo.crdTrue, pos);

				azel[1] = rtk->ssat[i].azel[1];
				azel[0] = rtk->ssat[i].azel[0];
				rp = re / (re + Hion) * cos(azel[1]);   // IPPz=asin(rp)
				ap = PI / 2.0 - azel[1] - asin(rp);
				sinap = sin(ap);
				tanap = tan(ap);
				cosaz = cos(azel[0]);
				posp[0] = asin(sin(pos[0]) * cos(ap) + cos(pos[0]) * sinap * cosaz);

				if ((pos[0] > 70.0 * D2R && tanap * cosaz > tan(PI / 2.0 - pos[0])) ||
					(pos[0]<-70.0 * D2R && -tanap * cosaz>tan(PI / 2.0 + pos[0]))) {
					posp[1] = pos[1] + PI - asin(sinap * sin(azel[0]) / cos(posp[0]));
				}
				else {
					posp[1] = pos[1] + asin(sinap * sin(azel[0]) / cos(posp[0]));
				}
			}
			else if (rtk->ssat[i].vsat[0] == 0) {
				ion = 99999.0;
				elev = 99999.0;
				vtec = 99999.0;
				posp[1] = 99999.0;
				posp[0] = 99999.0;
			}
		}
		else
		{
			ion = 99999.0;
			elev = 99999.0;
			vtec = 99999.0;
			posp[1] = 99999.0;
			posp[0] = 99999.0;
		}
		if (posp[1] == 99999.0) p += sprintf(p, "%9.3f%s%15.4f%s%15.4f%s", vtec, sep, posp[1], sep, posp[0], sep);
		else p += sprintf(p, "%9.3f%s%15.4f%s%15.4f%s", vtec, sep, posp[1] / D2R, sep, posp[0] / D2R, sep);
	}

	p += sprintf(p, "\n");

	n = p - (char*)buff;

	if (n > 0) {
		fwrite(buff, n, 1, fp);
	}
}
//output satellite elevation
static void outElev(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,n,i0,i1;
	char *p=(char *)buff;
	char *sep=" ";
	double sow,ep[6];
	int week;
	double elev=99999.0;

	time2epoch(time,ep);

	sow=time2gpst(time,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",(int)ep[0],sep,(int)ep[1],
		sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		if (rtk->ssat[i].vsat[0]==1) {
			elev=rtk->ssat[i].azel[1]*R2D;
		}
		else if (rtk->ssat[i].vsat[0]==0) {
			elev=99999.0;
		}
		p+=sprintf(p,"%9.3f%s",elev,sep);
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

//output geometry-free (GF) ambiguity
void outAmb_GF(FILE* fp, rtk_t *rtk, gtime_t t)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,n,i0,i1;
	char *p=(char *)buff;
    char *sep=" ";
	double sow,ep[6];
	int week;
	double gfm=0.0;

	time2epoch(t,ep);

	sow = time2gpst(t,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",
		(int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		if (rtk->ssat[i].vsat[0]==1) {
			gfm=rtk->ssat[i].gf;
		}
		else if (rtk->ssat[i].vsat[0]==0) {
			gfm=99999.0;
		}
		if (rtk->ssat[i].azel[1]<rtk->opt.elmin) {
			gfm=99999.0;
			p+=sprintf(p,"%9.3f%s",gfm,sep);
		}
		else
			if (gfm>10000.0) p+=sprintf(p,"%9.3f%s",gfm,sep);
			else p+=sprintf(p,"%9.3f%s",fmod(gfm,1000),sep);
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

/* output wide-lane ambiguity --------------------------------------------------
* input:  int i (0:before smoothing  1:after smoothing)
*-----------------------------------------------------------------------------*/
void outAmb_MW(FILE* fp, rtk_t *rtk, gtime_t t, int j)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,n,i0,i1;
	char *p=(char *)buff;
	char *sep=" ";
	double sow,ep[6];
	int week;
	double mwm=0.0;

	time2epoch(t,ep);

	sow = time2gpst(t,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",
		(int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		if (rtk->ssat[i].vsat[0]==1) {
			mwm=PPP_Glo.ssat_Ex[i].mw[j];
		}
		else if (rtk->ssat[i].vsat[0]==0) {
			mwm=99999.0;
		}
		if (rtk->ssat[i].azel[1]<rtk->opt.elmin) {
			mwm=99999.0;
			p+=sprintf(p,"%9.3f%s",mwm,sep);
		}
		else
			if (mwm>10000.0) p+=sprintf(p,"%9.3f%s",mwm,sep);
			else p+=sprintf(p,"%9.3f%s",fmod(mwm,1000),sep);
	}

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

/* output ionospheric-free (IF) ambiguity information for each satellite -------------------*/
void outAmb_IF(FILE* fp, rtk_t *rtk, gtime_t t)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,j=0,n,i0,i1;
	char *p=(char *)buff;
    char *sep=" ";
	double sow,ep[6];
	int week;
	double amb=0.0;

	time2epoch(t,ep);

	sow = time2gpst(t,&week);

	p+=sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",
		(int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,(int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,sow,sep);

	i0=MAXSAT; i1=0;
	if ((SYS_GPS&rtk->opt.navsys)) {
		if (0 < i0) i0=0;
		if (MAXPRNGPS > i1) i1=MAXPRNGPS;
	}
	if ((SYS_GLO&rtk->opt.navsys)) {
		if (MAXPRNGPS < i0) i0=MAXPRNGPS;
		if (MAXPRNGPS+MAXPRNGLO > i1) i1=MAXPRNGPS+MAXPRNGLO;
	}
	if ((SYS_GAL&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO < i0) i0=MAXPRNGPS+MAXPRNGLO;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
	}
	if ((SYS_QZS&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
	}
	if ((SYS_CMP&rtk->opt.navsys)) {
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS < i0) i0=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS;
		if (MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP > i1) i1=MAXPRNGPS+MAXPRNGLO+MAXPRNGAL+NSATQZS+MAXPRNCMP;
	}
	if (i0 > i1) { i0=0; i1=0; }

	for (i=i0;i<i1;i++) {
		j=IB(i+1,0,&rtk->opt);
		if (rtk->ssat[i].vsat[0]==1) {
			amb=rtk->x[j];//*FREQ1/CLIGHT;
		}
		else if (rtk->ssat[i].vsat[0]==0) {
			amb=99999.0;
		}
		p+=sprintf(p,"%9.3f%s",amb,sep);
    }

	p+=sprintf(p,"\n");

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}
}

//output initialized files for PPP in post-processing mode
static void outIppp(FILE *fp, rtk_t *rtk, gtime_t time)
{
	unsigned char buff[MAXSOLMSG+1];
	int i,n=0,sat;
	char *p=(char *)buff;
	double ep[6];
	char id[32];

	time2epoch(time,ep);

	for (i=0;i<MAXSAT;i++) {
		if (PPP_Info.ssat[i].vs==1) n++;
	}
	p+=sprintf(p,"> %4d %2d %2d %2d %2d %11.7f   %3d",(int)ep[0],
		(int)ep[1],(int)ep[2],(int)ep[3],(int)ep[4],ep[5],n);
	p+=sprintf(p,"\n");
	for (i=0;i<MAXSAT;i++) {
		if (PPP_Info.ssat[i].vs==0) continue;
		sat = i+1;
		satno2id(sat,id);
		p+=sprintf(p," %s,%d:%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,"
			"%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f,%15.4f",id,PPP_Info.ssat[sat-1].flag,
			PPP_Info.ssat[sat-1].satpos[0],PPP_Info.ssat[sat-1].satpos[1],PPP_Info.ssat[sat-1].satpos[2],
			PPP_Info.ssat[sat-1].satclk,PPP_Info.ssat[sat-1].azel[1]*R2D,PPP_Info.ssat[sat-1].azel[0]*R2D,
			PPP_Info.ssat[sat-1].P[0],PPP_Info.ssat[sat-1].P[1],PPP_Info.ssat[sat-1].P[2],
			PPP_Info.ssat[sat-1].L[0],PPP_Info.ssat[sat-1].L[1],PPP_Info.ssat[sat-1].L[2],
			PPP_Info.ssat[sat-1].shd,PPP_Info.ssat[sat-1].wmap,PPP_Info.ssat[sat-1].dsag,
			PPP_Info.ssat[sat-1].dtid,PPP_Info.ssat[sat-1].dant[0],PPP_Info.ssat[sat-1].dant[1],
			PPP_Info.ssat[sat-1].phw);
		p+=sprintf(p,"\n");
	}

	n=p-(char *)buff;

	if (n>0) {
		fwrite(buff,n,1,fp);
	}

}

//output result files
extern void outResult(rtk_t *rtk,const solopt_t *sopt)
{
	if (rtk->opt.ionoopt==IONOOPT_IF12||rtk->opt.ionoopt==IONOOPT_UC1) {
		if (PPP_Glo.outFp[OFILE_RESIC1]) outResi(PPP_Glo.outFp[OFILE_RESIC1],rtk,PPP_Glo.tNow,1,0,0);
		if (PPP_Glo.outFp[OFILE_RESIP1]) outResi(PPP_Glo.outFp[OFILE_RESIP1],rtk,PPP_Glo.tNow,1,1,0);
	}
	else if (rtk->opt.ionoopt==IONOOPT_UC12) {
		if (PPP_Glo.outFp[OFILE_RESIC1]) outResi(PPP_Glo.outFp[OFILE_RESIC1],rtk,PPP_Glo.tNow,1,0,0);
		if (PPP_Glo.outFp[OFILE_RESIP1]) outResi(PPP_Glo.outFp[OFILE_RESIP1],rtk,PPP_Glo.tNow,1,1,0);
		if (PPP_Glo.outFp[OFILE_RESIC2]) outResi(PPP_Glo.outFp[OFILE_RESIC2],rtk,PPP_Glo.tNow,2,0,0);
		if (PPP_Glo.outFp[OFILE_RESIP2]) outResi(PPP_Glo.outFp[OFILE_RESIP2],rtk,PPP_Glo.tNow,2,1,0);
	}
	if (PPP_Glo.outFp[OFILE_ELEV])  outElev(PPP_Glo.outFp[OFILE_ELEV],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_PDOP])  outPdop(PPP_Glo.outFp[OFILE_PDOP],rtk,PPP_Glo.tNow);
    if (PPP_Glo.outFp[OFILE_DTRP])  outDtrp(PPP_Glo.outFp[OFILE_DTRP],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_ISB])   outIsb(PPP_Glo.outFp[OFILE_ISB],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_ISB_M]) outIsb_m(PPP_Glo.outFp[OFILE_ISB_M],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_IFB])   outIfb(PPP_Glo.outFp[OFILE_IFB],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_LCPC])  outLCmPC(PPP_Glo.outFp[OFILE_LCPC],rtk,PPP_Glo.tNow);
	if (rtk->opt.ionoopt == IONOOPT_UC1 || rtk->opt.ionoopt == IONOOPT_UC12) {
		if (PPP_Glo.outFp[OFILE_STEC]) outStec(PPP_Glo.outFp[OFILE_STEC], rtk, PPP_Glo.tNow);
		if (PPP_Glo.outFp[OFILE_VTEC]) outVtec(PPP_Glo.outFp[OFILE_VTEC], rtk, PPP_Glo.tNow);
	}
	if (PPP_Glo.outFp[OFILE_AMBIF]) outAmb_IF(PPP_Glo.outFp[OFILE_AMBIF],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_AMBMW0]) outAmb_MW(PPP_Glo.outFp[OFILE_AMBMW0],rtk,PPP_Glo.tNow,0);
	if (PPP_Glo.outFp[OFILE_AMBMW1]) outAmb_MW(PPP_Glo.outFp[OFILE_AMBMW1],rtk,PPP_Glo.tNow,1);
	if (PPP_Glo.outFp[OFILE_AMBGF]) outAmb_GF(PPP_Glo.outFp[OFILE_AMBGF],rtk,PPP_Glo.tNow);
	if (PPP_Glo.outFp[OFILE_IPPP]) outIppp(PPP_Glo.outFp[OFILE_IPPP],rtk,PPP_Glo.tNow);
}
