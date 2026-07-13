Real-time STEC configuration files
==================================

Files copied here:
- igs20_2350.atx: satellite antenna phase-center model. Also usable for receiver antenna PCV if your antenna type exists in the file.
- ocnload_241227.blq: optional ocean loading BLQ example. For a user-owned receiver, use it only if it contains your station; otherwise generate a BLQ for your receiver position.
- rtstec_igs_rts.conf: real-time PPP/STEC template using receiver RTCM3 observations and IGS RTS/SSR RTCM3 corrections.

Fields you must edit before running:
- strpath[0]: your receiver observation stream.
- strpath[2]: your IGS RTS/SSR correction stream.
- sta_name: your receiver/station name.
- prcopt.anttype[0] and prcopt.posopt[1]: enable receiver antenna PCV only if the antenna type is known in the ATX file.

Output:
- out/*.pos: real-time PPP position solution.
- out/*.stat: real-time state output. With prcopt.ionoopt=4, STEC/slant ionosphere states are written as $ION records.

Notes:
- This repository's post-processing path writes .stec/.vtec files directly. The real-time server path writes ionosphere states through the status stream instead.
- prcopt.sateph=3 assumes APC-based IGS RTS/SSR corrections. If your SSR stream is COM-based, set prcopt.sateph=4.
- prcopt.tidecorr=1 is used by default because a user receiver usually has no fixed BLQ station record. Use tidecorr=7 only after setting filopt.blq to a BLQ file containing your station.
