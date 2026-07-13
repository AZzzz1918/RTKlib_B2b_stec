/*------------------------------------------------------------------------------
* results.c : PPP auxiliary result output functions
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define TECU_FACTOR(freq) (40.30E16/(freq)/(freq))
#define TEC_MISSING       99999.0
#define MIN_TEC_EL        (15.0*D2R)

/* number and index of PPP states, kept consistent with ppp.c ----------------*/
#define NF_(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP_(opt)     ((opt)->dynamics?9:3)
#ifdef BDS2BDS3
#define NC_(opt)     (NSYS+1)
#else
#define NC_(opt)     (NSYS)
#endif
#define ND_(opt)     (((opt)->nf>=3||(opt)->ionoopt==IONOOPT_EST)?NC_(opt):0)
#define NT_(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define II_(s,opt)   (NP_(opt)+NC_(opt)+ND_(opt)+NT_(opt)+(s)-1)

extern double Hion,re; /* set by iontec()/ionex.c, in km */

/* output satellite index range used by GAMP result files --------------------*/
static void outsatrange(int navsys, int *i0, int *i1)
{
    int gps0=0,gps1=NSATGPS;
    int glo0=gps1,glo1=glo0+NSATGLO;
    int gal0=glo1,gal1=gal0+NSATGAL;
    int qzs0=gal1,qzs1=qzs0+NSATQZS;
    int cmp0=qzs1,cmp1=cmp0+NSATCMP;

    *i0=MAXSAT;
    *i1=0;
    if (navsys&SYS_GPS) {if (gps0<*i0) *i0=gps0; if (gps1>*i1) *i1=gps1;}
    if (navsys&SYS_GLO) {if (glo0<*i0) *i0=glo0; if (glo1>*i1) *i1=glo1;}
    if (navsys&SYS_GAL) {if (gal0<*i0) *i0=gal0; if (gal1>*i1) *i1=gal1;}
    if (navsys&SYS_QZS) {if (qzs0<*i0) *i0=qzs0; if (qzs1>*i1) *i1=qzs1;}
    if (navsys&SYS_CMP) {if (cmp0<*i0) *i0=cmp0; if (cmp1>*i1) *i1=cmp1;}
    if (*i0>*i1) *i0=*i1=0;
}
/* write common time header --------------------------------------------------*/
static char *outtimehdr(char *p, gtime_t time, const char *sep)
{
    double sow,ep[6];
    int week;

    time2epoch(time,ep);
    sow=time2gpst(time,&week);
    sow=floor(sow*100.0+0.5)/100.0;
    if (sow>=604800.0) {
        sow-=604800.0;
        week++;
    }
    time2epoch(gpst2time(week,sow),ep);

    return p+sprintf(p,"%04d%s%02d%s%02d%s%02d%s%02d%s%02d%s%4d%s%9.2f%s",
                     (int)ep[0],sep,(int)ep[1],sep,(int)ep[2],sep,
                     (int)ep[3],sep,(int)ep[4],sep,(int)ep[5],sep,week,sep,
                     sow,sep);
}
/* L1/B1 factor for converting ionosphere delay to TECU ----------------------*/
static double tecfact(int sat, const ssat_t *ssat, const nav_t *nav)
{
    int prn=0,sys=ssat?satsys(sat,&prn):SYS_NONE;
    double freq=0.0;
    int i,fcn=0;

    if (ssat&&ssat->sys) sys=ssat->sys;
    if (sys==SYS_GLO) {
        if (nav) {
            if (0<prn&&prn<=32&&nav->glo_fcn[prn-1]>0) {
                fcn=nav->glo_fcn[prn-1]-8;
            }
            else {
                for (i=0;i<nav->ng;i++) {
                    if (nav->geph[i].sat==sat) {
                        fcn=nav->geph[i].frq;
                        break;
                    }
                }
            }
        }
        freq=FREQ1_GLO+DFRQ1_GLO*fcn;
    }
    else if (sys==SYS_CMP) {
        freq=FREQ1_CMP;
    }
    else {
        freq=FREQL1;
    }
    return TECU_FACTOR(freq);
}
/* GAMP-style L1/B1 factor for converting ionosphere delay to TECU -----------*/
static double gamp_tecfact(int sat, const ssat_t *ssat)
{
    int prn=0,sys=satsys(sat,&prn);
    double freq=0.0;

    if (ssat&&ssat->sys) sys=ssat->sys;
    if (sys==SYS_GLO) {
        freq=sat2freq(sat,CODE_L1C,NULL);
        if (freq<=0.0) freq=FREQ1_GLO;
    }
    else if (sys==SYS_CMP) {
        freq=FREQ1_CMP;
    }
    else {
        freq=FREQL1;
    }
    return TECU_FACTOR(freq);
}
/* slant ionosphere delay state in meters ------------------------------------*/
static int ionostate(const rtk_t *rtk, int sat, double *ion)
{
    int j;

    *ion=0.0;
    if (!rtk||rtk->opt.ionoopt!=IONOOPT_EST) return 0;
    j=II_(sat,&rtk->opt);
    if (j<0||j>=rtk->nx||rtk->x[j]==0.0) return 0;
    *ion=rtk->x[j];
    return 1;
}
/* output sTEC to buffer ------------------------------------------------------
* notes : PPP ionosphere states follow GAMP: slant L1 ionosphere delay.
*-----------------------------------------------------------------------------*/
extern int outstecs(uint8_t *buff, const rtk_t *rtk, const nav_t *nav,
                    gtime_t time)
{
    int i,i0,i1,j,n;
    char *p=(char *)buff;
    const char *sep=" ";
    double deg,stec,fact;

    if (!buff||!rtk) return 0;

    p=outtimehdr(p,time,sep);
    outsatrange(rtk->opt.navsys,&i0,&i1);

    for (i=i0;i<i1;i++) {
        stec=TEC_MISSING;
        j=II_(i+1,&rtk->opt);
        deg=rtk->ssat[i].azel[1]*R2D;
        if (deg>=15.0&&rtk->ssat[i].vsat[0]&&j>=0&&j<rtk->nx&&rtk->x[j]!=0.0) {
            fact=gamp_tecfact(i+1,rtk->ssat+i);
            if (fact>0.0) stec=fabs(rtk->x[j]/fact);
        }
        p+=sprintf(p,"%9.3f%s",stec,sep);
    }
    p+=sprintf(p,"\n");
    n=(int)(p-(char *)buff);
    return n>MAXSOLMSG?MAXSOLMSG:n;
}
/* output vTEC and IPP to buffer ---------------------------------------------
* notes : each satellite contributes vTEC, IPP longitude and IPP latitude.
*-----------------------------------------------------------------------------*/
extern int outvtecs(uint8_t *buff, const rtk_t *rtk, const nav_t *nav,
                    gtime_t time)
{
    int i,i0,i1,j,n;
    char *p=(char *)buff;
    const char *sep=" ";
    double pos[3]={0},posp[3]={TEC_MISSING,TEC_MISSING,0.0};
    double stec,vtec,fact,hion=Hion>0.0?Hion:HION/1000.0;
    double radius=re>0.0?re:RE_WGS84/1000.0;
    double deg,elev,mf,rp,ap,sinap,tanap,cosaz;
    const double hopt=506.7,a=0.9782,r=6378.137;

    if (!buff||!rtk) return 0;

    p=outtimehdr(p,time,sep);
    outsatrange(rtk->opt.navsys,&i0,&i1);
    ecef2pos(PPP_Glo.rr[0]!=0.0?PPP_Glo.rr:rtk->sol.rr,pos);

    for (i=i0;i<i1;i++) {
        vtec=TEC_MISSING;
        posp[0]=posp[1]=TEC_MISSING;
        j=II_(i+1,&rtk->opt);
        deg=rtk->ssat[i].azel[1]*R2D;
        if (deg>=15.0&&rtk->ssat[i].vsat[0]&&j>=0&&j<rtk->nx&&rtk->x[j]!=0.0) {
            fact=gamp_tecfact(i+1,rtk->ssat+i);
            if (fact>0.0) {
                stec=fabs(rtk->x[j]/fact);
                elev=PI/2.0-rtk->ssat[i].azel[1];
                mf=1.0/cos(asin(r*sin(a*elev)/(r+hopt)));
                vtec=stec/mf;
            }
            rp=radius/(radius+hion)*cos(rtk->ssat[i].azel[1]);
            ap=PI/2.0-rtk->ssat[i].azel[1]-asin(rp);
            sinap=sin(ap);
            tanap=tan(ap);
            cosaz=cos(rtk->ssat[i].azel[0]);
            posp[0]=asin(sin(pos[0])*cos(ap)+cos(pos[0])*sinap*cosaz);
            if ((pos[0]>70.0*D2R&&tanap*cosaz>tan(PI/2.0-pos[0]))||
                (pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI/2.0+pos[0]))) {
                posp[1]=pos[1]+PI-asin(sinap*sin(rtk->ssat[i].azel[0])/
                                        cos(posp[0]));
            }
            else {
                posp[1]=pos[1]+asin(sinap*sin(rtk->ssat[i].azel[0])/
                                     cos(posp[0]));
            }
        }
        if (posp[1]==TEC_MISSING) {
            p+=sprintf(p,"%9.3f%s%15.4f%s%15.4f%s",vtec,sep,posp[1],sep,
                       posp[0],sep);
        }
        else {
            p+=sprintf(p,"%9.3f%s%15.4f%s%15.4f%s",vtec,sep,posp[1]*R2D,sep,
                       posp[0]*R2D,sep);
        }
    }
    p+=sprintf(p,"\n");
    n=(int)(p-(char *)buff);
    return n>MAXSOLMSG?MAXSOLMSG:n;
}
/* output sTEC to file -------------------------------------------------------*/
extern void outstec(FILE *fp, const rtk_t *rtk, const nav_t *nav, gtime_t time)
{
    uint8_t buff[MAXSOLMSG+1];
    int n;

    if (!fp) return;
    if ((n=outstecs(buff,rtk,nav,time))>0) fwrite(buff,n,1,fp);
}
/* output vTEC to file -------------------------------------------------------*/
extern void outvtec(FILE *fp, const rtk_t *rtk, const nav_t *nav, gtime_t time)
{
    uint8_t buff[MAXSOLMSG+1];
    int n;

    if (!fp) return;
    if ((n=outvtecs(buff,rtk,nav,time))>0) fwrite(buff,n,1,fp);
}
