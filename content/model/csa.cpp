[PROB]
# Yoshikado et al. (2016)

- Title: __Quantitative Analyses of Hepatic OATP-Mediated
Interactions Between Statins and Inhibitors Using PBPK
Modeling With a Parameter Optimiaztion Method__
- Reference: CP\&T vol. 100 no. 5 pp. 513-23 11/2016
- Parameters: 40
- Compartments: 31

[CMT] 

// Dosing
igut

// CsA compartments
icent
me se ae
mc sc ac
iliv1 iliv2 iliv3 iliv4 iliv5

[PARAM] // CSA
iKp_mus = 2.98
iKp_adi = 17.3
iKp_ski = 13.6
iKp_liv = 16.7
ifb = 0.06


PSmus = 245/60
PSski = 37.4/60
PSadi = 10.2/60
ifhCLint = 0.587/60

ifafg = 0.572
iClr = 0
ika = 0.999
itlag =  0.254 

Vcent = 0.075


[PARAM]
Qh   = 1.200
Qmus = 0.642
Qski = 0.257
Qadi = 0.223

Vliv = 0.0241
Vmus = 0.4290
Vski = 0.1110
Vadi = 0.1430

exFliv = 0.278
exFmus = 0.146
exFski = 0.321
exFadi = 0.145


[MAIN]

if(NEWIND <=1) {
  double Vme = Vmus*exFmus;
  double Vae = Vadi*exFadi;
  double Vse = Vski*exFski;
  double Vmc = Vmus-Vme;
  double Vac = Vadi-Vae;
  double Vsc = Vski-Vse;
  double dVliv = Vliv/5.0;
}


[ODE]

// Statin concentrations

// CsA concentrations
double iCcent = icent/Vcent;

double Cme = me/Vme;
double Cse = se/Vse;
double Cae = ae/Vae;

double Cmc = mc/Vmc;
double Csc = sc/Vsc;
double Cac = ac/Vac;

double iCliv1 = iliv1/dVliv;
double iCliv2 = iliv2/dVliv;
double iCliv3 = iliv3/dVliv;
double iCliv4 = iliv4/dVliv;
double iCliv5 = iliv5/dVliv;

dxdt_igut = -ika/ifafg*igut;

dxdt_icent = 
  Qh*iCliv5/iKp_liv 
  - Qh*iCcent 
  - iClr*iCcent 
  - Qmus*(iCcent-Cme) 
  - Qski*(iCcent-Cse) 
  - Qadi*(iCcent-Cae);
  
dxdt_me = Qmus*(iCcent-Cme) - PSmus*ifb*(Cme-Cmc/iKp_mus);
dxdt_se = Qski*(iCcent-Cse) - PSski*ifb*(Cse-Csc/iKp_ski);
dxdt_ae = Qadi*(iCcent-Cae) - PSadi*ifb*(Cae-Cac/iKp_adi);
  
dxdt_mc = PSmus*ifb*(Cme-Cmc/iKp_mus);
dxdt_sc = PSski*ifb*(Cse-Csc/iKp_ski);
dxdt_ac = PSadi*ifb*(Cae-Cac/iKp_adi);
  
dxdt_iliv1 = Qh*(iCcent-iCliv1/iKp_liv) - (ifhCLint/5.0)*iCliv1 + ika*igut;
dxdt_iliv2 = Qh*(iCliv1-iCliv2)/iKp_liv - (ifhCLint/5.0)*iCliv2;
dxdt_iliv3 = Qh*(iCliv2-iCliv3)/iKp_liv - (ifhCLint/5.0)*iCliv3;
dxdt_iliv4 = Qh*(iCliv3-iCliv4)/iKp_liv - (ifhCLint/5.0)*iCliv4;
dxdt_iliv5 = Qh*(iCliv4-iCliv5)/iKp_liv - (ifhCLint/5.0)*iCliv5;
  


[TABLE]

capture CSA = icent/Vcent;
capture CSAliv = iliv1/dVliv;

