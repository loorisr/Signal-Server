#ifndef _ITM_TYPES_HH_
#define _ITM_TYPES_HH_

/* ITM state structures, shared by itwom3.0.cc and itm_batch.cc. */

struct prop_type {
	float aref;
	float dist;
	float hg[2];
	float rch[2];
	float wn;
	float dh;
	float dhd;
	float ens;
	float encc;
	float cch;
	float cd;
	float gme;
	float zgndreal;
	float zgndimag;
	float he[2];
	float dl[2];
	float the[2];
	float tiw;
	float ght;
	float ghr;
	float rph;
	float hht;
	float hhr;
	float tgh;
	float tsgh;
	float thera;
	float thenr;
	int rpl;
	int kwx;
	int mdp;
	int ptx;
	int los;
};

struct propv_type {
	float sgc;
	int lvar;
	int mdvar;
	int klim;
};

struct propa_type {
	float dlsa;
	float dx;
	float ael;
	float ak1;
	float ak2;
	float aed;
	float emd;
	float aes;
	float ems;
	float dls[2];
	float dla;
	float tha;
};

#endif /* _ITM_TYPES_HH_ */
