% This code compare genes regulated by information from different
% databases. 

% Written: 2020-05-28 Linnea Österberg 
% Updated:
respirationDataset1={};
respirationDataset2={};
respirationDataset3={};

respirationDataset1.Gis1upp="ADH1,TDH1,COX5B,PGK1,CIT1,ENO1,PDC1,FBP1,ERR1,ERR2,PCK1,ZWF1,PDX1,PDC5,PDC6,ICL1,MLS1,TKL2,GPP1,ERR3,CIT3,YJL045W,ALD2,GND2,TIM11,NDE2";
respirationDataset1.Gis1upp=split(respirationDataset1.Gis1upp,',');
respirationDataset1.Gis1upp=unique(respirationDataset1.Gis1upp);
respirationDataset1.Gis1down="QCR7,COX6,ATP2,COX8,HXK2,CYT1,COX9,COR1,RIP1,QCR8,ATP5,COX7,MDH1,KGD1,ATP15,SDH2,ATP7,NDI1,COX13,QCR10,ATP3,NDE1,GPD2,SOL4,ATP19,GPD1,SDH1,COX12,ATP17,ATP16,ATP20,ATP14";
respirationDataset1.Gis1down=split(respirationDataset1.Gis1down,',');
respirationDataset1.Gis1down=unique(respirationDataset1.Gis1down);
respirationDataset1.Msn2upp="TDH1,COX5B,CIT1,GAL10,HXK1,ATP4,PDC1,GAL7,ERR1,ERR2,COX7,PDX1,MDH1,GLK1,COX11,ATP15,SDH2,QCR9,ATP7,NDI1,TKL2,SDH3,PGM2,COX23,GPP2,ERR3,CIT3,ALD4,ALD2,PYK2,LSC2,GPD1,HFD1,COX20,ATP20";
respirationDataset1.Msn2upp=split(respirationDataset1.Msn2upp,',');
respirationDataset1.Msn2upp=unique(respirationDataset1.Msn2upp);
respirationDataset1.Msn2down="TDH2,CDC19,PGK1,ENO1,GPM1,COX4,HXK2,COX9,COR1,PGI1,FBA1,PFK1,PFK2,PGM1,GND1,NDE1,ACS2,SOL4,FAP7";
respirationDataset1.Msn2down=split(respirationDataset1.Msn2down,',');
respirationDataset1.Msn2down=unique(respirationDataset1.Msn2down);
respirationDataset1.Msn4upp="COX5B,HXK1,PDC1,COX7,PDX1,MDH1,GLK1,COX11,SDH2,ATP7,NDI1,TKL2,SDH3,PGM2,COX23,GPP2,ERR3,CIT3,ALD4,ALD2,PYK2,LSC2,GPD1,HFD1,COX20,NDE2,ATP20";
respirationDataset1.Msn4upp=split(respirationDataset1.Msn4upp,',');
respirationDataset1.Msn4upp=unique(respirationDataset1.Msn4upp);
respirationDataset1.Msn4down="CDC19,GPM1,COX9,PGI1,GPP1,SOL4,RKI1";
respirationDataset1.Msn4down=split(respirationDataset1.Msn4down,',');
respirationDataset1.Msn4down=unique(respirationDataset1.Msn4down);
respirationDataset1.Adr1upp="ADH1,TDH2,TDH3,CDC19,PGK1,ENO2,ENO1,PDC1,FBP1,ERR1,COX7,PYC1,FBA1,PDC5,TKL2,ALD5,GPD2,ERR3,CIT3,ALD4,GND2,LSC1,ACS1,HFD1";
respirationDataset1.Adr1upp=split(respirationDataset1.Adr1upp,',');
respirationDataset1.Adr1upp=unique(respirationDataset1.Adr1upp);
respirationDataset1.Adr1down="COX5B,COX4,COX8,ATP5,QCR9,QCR10,ALD6,ATP19";
respirationDataset1.Adr1down=split(respirationDataset1.Adr1down,',');
respirationDataset1.Adr1down=unique(respirationDataset1.Adr1down);
respirationDataset1.Sip4upp="ACS1";
respirationDataset1.Sip4down="PYC1,PGI1,NDE1";
respirationDataset1.Sip4down=split(respirationDataset1.Sip4down,',');
respirationDataset1.Cat8upp="PGK1,ATP2,FUM1,FBP1,PCK1,ICL1,ALD4,YJL045W,ALD6,SDH1,ACS1,RKI1";
respirationDataset1.Cat8upp=split(respirationDataset1.Cat8upp,',');
respirationDataset1.Cat8upp=unique(respirationDataset1.Cat8upp);
respirationDataset1.Cat8down= "PGI1";


respirationDataset2.Gis1upp={};
respirationDataset2.Gis1down={};
respirationDataset2.Msn2upp="ADP1,AGE1,AIM17,ASK1,ASP3-1,ASP3-2,ASP3-3,ASP3-4,AVT1,BMH1,CCP1,CIP1,CMK2,COQ6,CTT1,CUE4,CUZ1,CWP1,DDR48,DIA1,EFB1,ENV11,FMP48,GIS4,GND2,GPH1,GPT2,GRE3,GSY2,HER1,HPF1,HSP30,HSP42,HSP78,ISU2,LST8,MDH1,MF(ALPHA)2,MGA1,MRH1,MRS4,MRX16,MSN4,NSG2,PDC6,PEX32,PGK1,POP7,PST2,PTH4,RBA50,ROD1,RRT15,SFK1,SIS1,SKM1,SNA2,SNR18,SSA1,SSC1,TAH11,TDH3,tE(UUC)E1,tG(GCC)P2,TOS8,TPI1,TPK1,TPO4,TSL1,TVP23,UPF3,YBR085C-A,YDR524C-B,YDR524W-C,YGP1,YKL063C,YKL096C-B,YLR159W,YLR162W,YML100W-A,YMR316C-A,YOR385W,YPR159C-A";
respirationDataset2.Msn2upp=split(respirationDataset2.Msn2upp,',');
respirationDataset2.Msn2upp=unique(respirationDataset2.Msn2upp);
respirationDataset2.Msn2down={};
respirationDataset2.Msn4upp="AFR1,ARF1,ASI1,CCP1,CTT1,CUE4,DDR48,EFB1,ENV11,FMP48,GPT2,GRE3,HER1,HOR7,HSP26,HSP30,HSP78,ISU2,LSM3,LST8,MCP1,MDH1,MMF1,MRPL4,MRS4,NNR1,PCL7,PDC6,PEX32,POP7,PRE2,ROD1,ROX1,RPL11A,RTC3,SIS1,SNR18,SSA1,SSA4,TDH3,TSA2,TSL1,TVP23,UPF3,UTP11,WTM2,YGP1,YHR033W,YHR086W-A,YML100W-A,YPR063C,ZWF1";
respirationDataset2.Msn4upp=split(respirationDataset2.Msn4upp,',');
respirationDataset2.Msn4upp=unique(respirationDataset2.Msn4upp);
respirationDataset2.Msn4down={};
respirationDataset2.Adr1upp="USV1,ATG32,ATP7,CMK2,CYK3,FUN26,GLC7,GLD1,HEM25,HIR3,IPA1,MIC27,NCW2,PPT2,PUT3,PXA1,SNR52,tT(UGU)P,VPS27,YPR108W-A";
respirationDataset2.Adr1upp=split(respirationDataset2.Adr1upp,',');
respirationDataset2.Adr1upp=unique(respirationDataset2.Adr1upp);
respirationDataset2.Adr1down={};
respirationDataset2.Sip4upp="CMK2,GCV1,ICL1,MDH2,OYE2,PCK1,PML39,RRT15,TRM13,YLR162W,FBP1";
respirationDataset2.Sip4upp=split(respirationDataset2.Sip4upp,',');
respirationDataset2.Sip4upp=unique(respirationDataset2.Sip4upp);
respirationDataset2.Sip4down={};
respirationDataset2.Cat8upp="FBP1,ICL1,PCK1,AIM17,FRE1,HSP104,HSP30,MET17,PHO89,TSL1,YKL177W";
respirationDataset2.Cat8upp=split(respirationDataset2.Cat8upp,',');
respirationDataset2.Cat8upp=unique(respirationDataset2.Cat8upp);
respirationDataset2.Cat8down={};


respirationDataset3.Gis1upp="FAA2,IDP2,PCK1,ICL1,FBP1,MLS1";
respirationDataset3.Gis1upp=split(respirationDataset3.Gis1upp,',');
respirationDataset3.Gis1upp=unique(respirationDataset3.Gis1upp);
respirationDataset3.Gis1down={};
respirationDataset3.Msn2upp="GDH3,TKL2,GLK1,ARA1,SSA3,HXK1,HSP104,HOR2,SSA4,PGM2,SOD2,TPS1,ALD3,GPD1,MDH2,PYK2,GUT2,CTA1,TDH3,GPM1,TPI1,FBA1,MDH1,GID8";
respirationDataset3.Msn2upp=split(respirationDataset3.Msn2upp,',');
respirationDataset3.Msn2upp=unique(respirationDataset3.Msn2upp);
respirationDataset3.Msn2down="ICL1,ADH2,ALD6,KGD2";
respirationDataset3.Msn2down=split(respirationDataset3.Msn2down,',');
respirationDataset3.Msn2down=unique(respirationDataset3.Msn2down);
respirationDataset3.Msn4upp="GDH3,TKL2,ARA1,SSA3,HSP104,HOR2,SSA4,PGM2,SOD2,TPS1,ALD3,GPD1,MDH2,PYK2,HXK1,GLK1,GUT2,ADR1,MDH1,TDH3";
respirationDataset3.Msn4upp=split(respirationDataset3.Msn4upp,',');
respirationDataset3.Msn4upp=unique(respirationDataset3.Msn4upp);
respirationDataset3.Msn4down="ICL1,ADH2,ALD6,KGD2";
respirationDataset3.Msn4down=split(respirationDataset3.Msn4down,',');
respirationDataset3.Msn4down=unique(respirationDataset3.Msn4down);
respirationDataset3.Adr1upp="FOX2,DCI1,POX1,SPS19,CTA1,CYB2,ADH2,GUT2,POT1,ALD4ICL2,CIT3,PXA1,GUT1,IDP3,PEX11,MDH2,MLS1,FBP1,ADY2,ADH2,ALD4,FDH1,POX1,FOX2,ATO3,ICL1,JEN1,CYB2,SPS19,SPS19,FDH1,FOX2,SPG1,CYB2,CTA1,FBP1,ADY2,POX1,ALD4,ADH2,ADH2,ADY2,CTA1,POT1,FOX2,POX1,ADH2,ADH2,ACS1,ADY2,ALD4,LSC2,POT1,CTA1,JEN1,FOX2,CIT3,GUT1,ICL2,GUT2,ADH2,CYB2,CAT2,IDP3,HXT6,ACO1,PGI1,PDC1,KGD2,STL1,JEN1,FBP1,ADH2,PEX25,CIT1,FOX2,CYB2,POX1,CTA1,GUT2,CIT3,POX1,ALD4,ACS1,CTA1,POT1,GUT1,ADH2,ICL2,CTA1,POX1,POT1,FOX2,SPS19,PIP2,FAA1,POX1,CTA1,PEX11,POT1,SPS19,FOX2,PXA1";
respirationDataset3.Adr1upp=split(respirationDataset3.Adr1upp,',');
respirationDataset3.Adr1upp=unique(respirationDataset3.Adr1upp);
respirationDataset3.Adr1down="POT1,CTA1,SPG1,MDH2";
respirationDataset3.Adr1down=split(respirationDataset3.Adr1down,',');
respirationDataset3.Adr1down=unique(respirationDataset3.Adr1down);
respirationDataset3.Sip4upp="FBP1,ICL1,PCK1,MDH2,SFC1,ICL1,PCK1,MLS1,MLS1,SFC1,PCK1,FBP1,MDH2,ICL1";
respirationDataset3.Sip4upp=split(respirationDataset3.Sip4upp,',');
respirationDataset3.Sip4upp=unique(respirationDataset3.Sip4upp);
respirationDataset3.Sip4down={};
respirationDataset3.Cat8upp="FBP1,ADH2,ICL1,MDH2,IDP2,MLS1,MLS1,ADY2,FBP1,SFC1,FBP1,CAT2,ACS1,YAT2,ADH2,MDH2,ICL1,MLS1,ALD6,DLD1,ACH1,JEN1,PCK1,REG2,IDP2,SFC1,ICL1,FBP1,PCK1,MDH2,MLS1,ACS1,MDH2,MLS1,FBP1,ADY2,ADH2,ACS1,ALD4,FDH1,POX1,FOX2,ATO3,ICL1,JEN1,CYB2,SPS19,JEN1,ATO3,POT1,ACS1,ICL2,IDP2,REG2,ICL1,GUT1,CTA1,POX1,JEN1,PCK1,FBP1,YAT2,SFC1,ADY2,POX1,ADH2,ATO3,MDH2,ACS1,MLS1,SPS19,JEN1,FDH1";
respirationDataset3.Cat8upp=split(respirationDataset3.Cat8upp,',');
respirationDataset3.Cat8upp=unique(respirationDataset3.Cat8upp);
respirationDataset3.Cat8down="POT1,CTA1,SPG1,POT1,ALD4,CYB2,CTA1,SPG1";
respirationDataset3.Cat8down=split(respirationDataset3.Cat8down,',');
respirationDataset3.Cat8down=unique(respirationDataset3.Cat8down);

fermenrtationDataset1={};
fermenrtationDataset2={};
fermenrtationDataset3={};

fermenrtationDataset1.Mig1upp="QCR6,ADH1,TDH3,TDH1,COX6,CDC19,PGK1,ENO2,ENO1,TPI1,FBA1,ALD5,SOL4,GND2,ALD6,GPD1,COX12";
fermenrtationDataset1.Mig1upp=split(fermenrtationDataset1.Mig1upp,',');
fermenrtationDataset1.Mig1upp=unique(fermenrtationDataset1.Mig1upp);
fermenrtationDataset1.Mig1down="GAL1,GAL10,ATP4,FBP1,TIM11";
fermenrtationDataset1.Mig1down=split(fermenrtationDataset1.Mig1down,',');
fermenrtationDataset1.Mig1down=unique(fermenrtationDataset1.Mig1down);
fermenrtationDataset1.Sfp1upp="ADH1,TDH2,TDH3,PGK1,ENO2,ENO1,GPM1,GAL1,GAL10,HXK2,PDC1,GAL7,ERR1,ZWF1,LAT1,PGI1,FBA1,TAL1,PDC5,PFK2,GLK1,ACO1,TKL1,PDC6,GND1,SOL3,ALD5,GPP1,ERR3,COX18,TIM11,HFD1,YLR446W,RKI1,COX19";
fermenrtationDataset1.Sfp1upp=split(fermenrtationDataset1.Sfp1upp,',');
fermenrtationDataset1.Sfp1upp=unique(fermenrtationDataset1.Sfp1upp);
fermenrtationDataset1.Sfp1down="QCR6,COX6,CDC19,CIT1,COX4,COX8,ATP4,COX9,COR1,QCR2,RIP1,YCR8,PCK1,MDH1,ATP15,SDH2,QCR9,ATP7,COX13,SDH3,QCR10,GPP2,ALD4,LSC2,ATP18,ATP17,NDE2,FAP7,ATP16";
fermenrtationDataset1.Sfp1down=split(fermenrtationDataset1.Sfp1down,',');
fermenrtationDataset1.Sfp1down=unique(fermenrtationDataset1.Sfp1down);

fermenrtationDataset2.Mig1upp={};
fermenrtationDataset2.Mig1upp=split(fermenrtationDataset2.Mig1upp,',');
fermenrtationDataset2.Mig1upp=unique(fermenrtationDataset2.Mig1upp);
fermenrtationDataset2.Mig1down="ARC1,ENA1,GAL1,GAL2,GAL4,HXT2,HXT3,HXT4,MAL61,MAL62,MAL63,NYV1,SNF3,SUC2,AGA1,AIM17,ALD6,APL1,APT2,BTN2,DDR48,ESP1,FRE1,GAC1,GDB1,HXT3,HXT6,MAL33,MBF1,MPC3,NCW2,ODC1,PHO89,PPH3,RDN25-1,REG2,RRT15,STP4,YDL023C,YDR444W,YGP1,YLR162W,YMR085W,YOL014W";
fermenrtationDataset2.Mig1down=split(fermenrtationDataset2.Mig1down,',');
fermenrtationDataset2.Mig1down=unique(fermenrtationDataset2.Mig1down);
fermenrtationDataset2.Sfp1upp="AAH1,ACM1,ACS2,ADE12,ADE17,ADE3,ADH1,ADH2,ADH3,ADH5,ADK1,AGA2,AGP2,AHP1,AIM17,ALD5,ALD6,ANB1,APE3,APE4,APS1,APT1,AQR1,ARA2,ARB1,ARC1,ARO4,ARO9,ASC1,ASE1,ASN2,ATF2,ATP20,ATP7,BAT1,BAT2,BNI5,BRR1,BUD28,CAB3,CAM1,CAR1,CAR2,CBF5,CCE1,CCL1,CDC5,CHM7,CIP1,CLB1,CLN3,CMK2,COX7,CSM3,CSN9,CSS1,CTO1,CTS1,CUB1,CYC1,CYS3,CYS4,CYT1,DBF4,DDR48,DED81,DFG10,DIA1,DIC1,DPC25,DPS1,ECM13,EEB1,EFB1,EFM1,EFT2,EGD1,ELM1,ELO3,ENO1,ENO2,EPL1,ERB1,ERC1,ERD2,ERR2,ERR3,ERV14,ERV2,ESA1,EXG1,FBA1,FCP1,FCY2,FDO1,FET3,FET4,FPR4,FRD1,FRS1,FTR1,FUN12,FUR1,FUS1,FUS2,GAR1,GAT2,GCD11,GCV1,GCV3,GDS1,GDT1,GGC1,GIC2,GIS2,GLG2,GLK1,GPA1,GPB2,GPD2,GPM1,GPP1,GSC2,GSP1,GUS1,GUT2,GYP8,HEF3,HIS3,HOL1,HOM2,HOM3,HOP2,HPF1,HPT1,HSH49,HSP12,HSP30,HSP60,HSP82,HTC1,HXK1,HXK2,HXT16,HYP2,HYR1,ICS2,ICS3,ILV2,ILV3,ILV5,ILV6,IMD2,IMD3,IMD4,IPI3,IZH4,JHD2,KAP95,KAR4,KAR9,KNS1,KRE29,LCB1,LEA1,LHP1,LIA1,LIP5,LOH1,LRP1,LSB1,MAE1,MAP1,MBF1,MCM16,MDH2,MET6,MFA1,MFM1,MIC26,MIG2,MMF1,MRM2,MRP51,MSN4,NCE103,NDC80,NDE1,NEW1,NFU1,NIP1,NIP100,NOG2,NOP1,NOP56,NOP58,NRE1,NRG2,NRP1,NSR1,NTC20,NTF2,NTG2,NTR2,OKP1,OLA1,OPI11,ORC3,PCS60,PDC1,PDC5,PDC6,PDE2,PDI1,PDR8,PER33,PET20,PEX27,PEX30,PGI1,PGK1,PHO11,PHO12,PHO3,PHO5,PHO8,PHO89,PLB2,PMT1,POL4,PPT1,PRE2,PRM1,PRP40,PRP43,PRS1,PRS3,PRS5,PRT1,PRY1,PSE1,PUL4,RBG1,RCL1,RDL1,RFM1,RFU1,RGD2,RIE1,RIX7,RKI1,RLI1,RPA135,RPB5,RPB8,RPL10,RPL11A,RPL11B,RPL12A,RPL12B,RPL13A,RPL13B,RPL14A,RPL14B,RPL15A,RPL15B,RPL16A,RPL16B,RPL17A,RPL17B,RPL18A,RPL18B,RPL1A,RPL1B,RPL20A,RPL20B,RPL21A,RPL21B,RPL22A,RPL22B,RPL23A,RPL23B,RPL24A,RPL25,RPL26A,RPL26B,RPL27A,RPL28,RPL29,RPL2A,RPL2B,RPL3,RPL30,RPL31A,RPL32,RPL33A,RPL33B,RPL34A,RPL34B,RPL36A,RPL36B,RPL37A,RPL38,RPL40A,RPL42A,RPL42B,RPL43A,RPL43B,RPL4A,RPL4B,RPL5,RPL6A,RPL6B,RPL7A,RPL7B,RPL8A,RPL8B,RPL9B,RPM2,RPP0,RPP1A,RPP1B,RPP2A,RPS0A,RPS0B,RPS10A,RPS10B,RPS11A,RPS11B,RPS12,RPS13,RPS14A,RPS15,RPS16B,RPS17A,RPS17B,RPS18A,RPS18B,RPS19A,RPS19B,RPS1A,RPS1B,RPS2,RPS20,RPS21A,RPS21B,RPS22A,RPS22B,RPS23A,RPS24B,RPS25A,RPS25B,RPS26A,RPS26B,RPS28B,RPS29A,RPS29B,RPS3,RPS30A,RPS30B,RPS31,RPS4A,RPS4B,RPS5,RPS6A,RPS6B,RPS7B,RPS8A,RPS8B,RPS9A,RPS9B,RRN11,RRT15,RRT2,RTN1,RXT2,SCP160,SCW11,SCW4,SED1,SES1,SFI1,SHM2,SIR1,SKS1,SLI15,SLK19,SMC6,SND2,SND3,SNR17A,SNR44,SNR58,SNR59,SNU13,SNZ1,SPE3,SPP1,SRL2,SSA2,SSB1,SSB2,SSC1,SSU1,STB6,STE2,STI1,STM1,SUI2,SUP45,SYO1,TAM41,TDH1,TDH2,TEF4,TFG2,THR1,THR4,TIF1,TIF34,TIF5,TIF6,TIM23,TKL1,TMA19,TMA46,TPM1,TRM112,TRX1,TSA1,TYS1,ULP1,ULS1,UPS2,URA5,UTP15,UTP4,VAS1,VID24,VPH1,VPS20,VPS24,VPS45,WRS1,WSC4,YAH1,YAL004W,YAR075W,YBL077W,YBR063C,YBR137W,YBR287W,YCR013C,YCR087C-A,YCR102C,YDL228C,YDL241W,YDR154C,YDR341C,YDR417C,YEF3,YER158C,YFL064C,YGL102C,YGP1,YGR137W,YGR160W,YHB1,YHL009W-A,YHL049C,YHR020W,YIL082W,YJL114W,YKL030W,YKL123W,YKL153W,YLL066C,YLR041W,YLR076C,YLR156W,YLR162W,YLR171W,YLR225C,YLR339C,YLR349W,YLR446W,YLR460C,YLR462W,YML6,YMR075C-A,YMR085W,YMR119W-A,YMR173W-A,YNL146W,YOL014W,YOP1,YOR296W,YOR338W,YOR385W,YPL197C,YPL245W,YPL251W,YPR003C,YPS1,YPT52,YRA1,YRF1-5,ZEO1,ZPR1,ZRT1,ZWF1";
fermenrtationDataset2.Sfp1upp=split(fermenrtationDataset2.Sfp1upp,',');
fermenrtationDataset2.Sfp1upp=unique(fermenrtationDataset2.Sfp1upp);
fermenrtationDataset2.Sfp1down={};


fermenrtationDataset3.Mig1upp={};%"CAT8,FBP1,HXT1,HXT3,PCK1,HXT2,CAT8,HXK1,GLK1,FBP1,HXT4,SNF3,HAP4,PDC1,MLS1,ICL1,MDH2,HXT1,HXT3";
fermenrtationDataset3.Mig1upp=split(fermenrtationDataset3.Mig1upp,',');
fermenrtationDataset3.Mig1upp=unique(fermenrtationDataset3.Mig1upp);
fermenrtationDataset3.Mig1down="MDH2,JEN1,ICL1,HXT3,CRC1,HXK1,LSC2,GUT1,HXT1,GUT2,CAT8,ADR1,HXT4,REG2,HXT2,MDH3,SFC1,GUT2,MLS1,SNF3,HXT4,MIG1,HAP4,HXT2,ICL1,HXT3,HXK1,PDC1,FBP1,CAT8,GLK1,PCK1,PCK1,FBP1,PDC1,HAP4,FBP1,CAT8,GUT1,GUT2,ADR1,HXK1,JEN1,REG2,GUT1,LSC2,HXT4,CRC1,ICL1,CAT8,MDH2,HXT3,HXT2,HXT2,HXT6,HXT4,HXT7,ICL1,REG2,CAT8,ADR1,CRC1,JEN1,GUT1,HXT2,LSC2,GUT2,HXK1,HXT4";
fermenrtationDataset3.Mig1down=split(fermenrtationDataset3.Mig1down,',');
fermenrtationDataset3.Mig1down=unique(fermenrtationDataset3.Mig1down);
fermenrtationDataset3.Sfp1upp="ADH2";
fermenrtationDataset3.Sfp1upp=split(fermenrtationDataset3.Sfp1upp,',');
fermenrtationDataset3.Sfp1upp=unique(fermenrtationDataset3.Sfp1upp);
fermenrtationDataset3.Sfp1down={};

% Two TF are according to the same source both upregulating and
% downregulating the same gene. Probably inpropperly annotated about the
% conditions. 
C = intersect(respirationDataset3.Adr1upp,respirationDataset3.Adr1down);
disp('According to dataset 3, these genes are both uppregulated and downregulated by Adr1:');
disp( C);
C = intersect(respirationDataset3.Cat8upp,respirationDataset3.Cat8down);
disp('According to dataset 3, these genes are both uppregulated and downregulated by Cat8:');
disp( C);

%% Intersect all data: 

%merge upregulated and downregulted genes for all datasets
respirationDataset1upp=union(respirationDataset1.Gis1upp, respirationDataset1.Msn2upp);
respirationDataset1upp=union(respirationDataset1upp, respirationDataset1.Msn4upp);
respirationDataset1upp=union(respirationDataset1upp, respirationDataset1.Adr1upp);
respirationDataset1upp=union(respirationDataset1upp, respirationDataset1.Sip4upp);
respirationDataset1upp=union(respirationDataset1upp, respirationDataset1.Cat8upp);

respirationDataset1down=union(respirationDataset1.Gis1down, respirationDataset1.Msn2down);
respirationDataset1down=union(respirationDataset1down, respirationDataset1.Msn4down);
respirationDataset1down=union(respirationDataset1down, respirationDataset1.Adr1down);
respirationDataset1down=union(respirationDataset1down, respirationDataset1.Sip4down);
respirationDataset1down=union(respirationDataset1down, respirationDataset1.Cat8down);

respirationDataset2upp=union(respirationDataset2.Gis1upp, respirationDataset2.Msn2upp);
respirationDataset2upp=union(respirationDataset2upp, respirationDataset2.Msn4upp);
respirationDataset2upp=union(respirationDataset2upp, respirationDataset2.Adr1upp);
respirationDataset2upp=union(respirationDataset2upp, respirationDataset2.Sip4upp);
respirationDataset2upp=union(respirationDataset2upp, respirationDataset2.Cat8upp);

respirationDataset2down=union(respirationDataset2.Gis1down, respirationDataset2.Msn2down);
respirationDataset2down=union(respirationDataset2down, respirationDataset2.Msn4down);
respirationDataset2down=union(respirationDataset2down, respirationDataset2.Adr1down);
respirationDataset2down=union(respirationDataset2down, respirationDataset2.Sip4down);
respirationDataset2down=union(respirationDataset2down, respirationDataset2.Cat8down);

respirationDataset3upp=union(respirationDataset3.Gis1upp, respirationDataset3.Msn2upp);
respirationDataset3upp=union(respirationDataset3upp, respirationDataset3.Msn4upp);
respirationDataset3upp=union(respirationDataset3upp, respirationDataset3.Adr1upp);
respirationDataset3upp=union(respirationDataset3upp, respirationDataset3.Sip4upp);
respirationDataset3upp=union(respirationDataset3upp, respirationDataset3.Cat8upp);

respirationDataset3down=union(respirationDataset3.Gis1down, respirationDataset3.Msn2down);
respirationDataset3down=union(respirationDataset3down, respirationDataset3.Msn4down);
respirationDataset3down=union(respirationDataset3down, respirationDataset3.Adr1down);
respirationDataset3down=union(respirationDataset3down, respirationDataset3.Sip4down);
respirationDataset3down=union(respirationDataset3down, respirationDataset3.Cat8down);

fermentationDataset1upp=union(fermenrtationDataset1.Mig1upp,fermenrtationDataset1.Sfp1upp);
fermentationDataset2upp=union(fermenrtationDataset2.Mig1upp,fermenrtationDataset2.Sfp1upp);
fermentationDataset3upp=union(fermenrtationDataset3.Mig1upp,fermenrtationDataset3.Sfp1upp);

fermentationDataset1down=union(fermenrtationDataset1.Mig1down,fermenrtationDataset1.Sfp1down);
fermentationDataset2down=union(fermenrtationDataset2.Mig1down,fermenrtationDataset2.Sfp1down);
fermentationDataset3down=union(fermenrtationDataset3.Mig1down,fermenrtationDataset3.Sfp1down);

% Look at intersects 
respirationI12upp=intersect(respirationDataset1upp,respirationDataset2upp);
respirationI13upp=intersect(respirationDataset1upp,respirationDataset3upp);
respirationI23upp=intersect(respirationDataset2upp,respirationDataset3upp);
respirationI123upp=intersect(respirationI12upp,respirationI23upp);
respirationJoinupp=union(respirationI12upp, respirationI13upp);
respirationJoinupp=union(respirationJoinupp, respirationI23upp);

respirationI12down=intersect(respirationDataset1down,respirationDataset2down);
respirationI13down=intersect(respirationDataset1down,respirationDataset3down);
respirationI23down=intersect(respirationDataset2down,respirationDataset3down);
respirationI123down=intersect(respirationI12down,respirationI23down);
respirationJoindown=union(respirationI12down, respirationI13down);
respirationJoindown=union(respirationJoindown, respirationI23down);

fermentationI12upp=intersect(fermentationDataset1upp,fermentationDataset2upp);
fermentationI13upp=intersect(fermentationDataset1upp,fermentationDataset3upp);
fermentationI23upp=intersect(fermentationDataset2upp,fermentationDataset3upp);
fermentationI123upp=intersect(fermentationI12upp,fermentationI23upp);
fermentationJoinupp=union(fermentationI12upp, fermentationI13upp);
fermentationJoinupp=union(fermentationJoinupp, fermentationI23upp);

fermentationI12down=intersect(fermentationDataset1down,fermentationDataset2down);
fermentationI13down=intersect(fermentationDataset1down,fermentationDataset3down);
fermentationI23down=intersect(fermentationDataset2down,fermentationDataset3down);
fermentationI123down=intersect(fermentationI12down,fermentationI23down);
fermentationJoindown=union(fermentationI12down, fermentationI13down);
fermentationJoindown=union(fermentationJoindown, fermentationI23down);

%% Look at overlap in model
load('../../../models/reduced_ecYeast_fermentation.mat');
model=table();
model.enzymes      = ecModel_ferm.enzymes;
model.genes        = ecModel_ferm.enzNames;
pathways           = mapEnzymeSubSystems(ecModel_ferm.enzymes,ecModel_ferm);
model.subSystems   = pathways;
model.respirationDataset1upp = zeros(length(model.enzymes),1);
model.respirationDataset1down = zeros(length(model.enzymes),1);
model.respirationDataset2upp = zeros(length(model.enzymes),1);
model.respirationDataset2down = zeros(length(model.enzymes),1);
model.respirationDataset3upp = zeros(length(model.enzymes),1);
model.respirationDataset3down = zeros(length(model.enzymes),1);
model.respirationJoinupp = zeros(length(model.enzymes),1);
model.respirationJoindown = zeros(length(model.enzymes),1);
model.fermentationDataset1upp = zeros(length(model.enzymes),1);
model.fermentationDataset1down = zeros(length(model.enzymes),1);
model.fermentationDataset2upp = zeros(length(model.enzymes),1);
model.fermentationDataset2down = zeros(length(model.enzymes),1);
model.fermentationDataset3upp = zeros(length(model.enzymes),1);
model.fermentationDataset3down = zeros(length(model.enzymes),1);
model.fermentationJoinupp = zeros(length(model.enzymes),1);
model.fermentationJoindown = zeros(length(model.enzymes),1);
for i = 1:length(model.enzymes)
    model.respirationDataset1upp(i)= sum(ismember(respirationDataset1upp,model.genes(i)));
    model.respirationDataset1down(i) = sum(ismember(respirationDataset1down,model.genes(i)));
    model.respirationDataset2upp(i) = sum(ismember(respirationDataset2upp,model.genes(i)));
    model.respirationDataset2down(i) = sum(ismember(respirationDataset2down,model.genes(i)));
    model.respirationDataset3upp(i) = sum(ismember(respirationDataset3upp,model.genes(i)));
    model.respirationDataset3down(i) = sum(ismember(respirationDataset3down,model.genes(i)));
    model.fermentationDataset1upp(i) = sum(ismember(fermentationDataset1upp,model.genes(i)));
    model.fermentationDataset1down(i) = sum(ismember(fermentationDataset1down,model.genes(i)));
    model.fermentationDataset2upp(i) = sum(ismember(fermentationDataset2upp,model.genes(i)));
    model.fermentationDataset2down(i) = sum(ismember(fermentationDataset2down,model.genes(i)));
    model.fermentationDataset3upp(i) = sum(ismember(fermentationDataset3upp,model.genes(i)));
    model.fermentationDataset3down(i) = sum(ismember(fermentationDataset3down,model.genes(i)));
    model.respirationJoinupp(i) = sum(ismember(respirationJoinupp,model.genes(i)));
    model.respirationJoindown(i) = sum(ismember(respirationJoindown,model.genes(i)));
    model.fermentationJoinupp(i) = sum(ismember(fermentationJoinupp,model.genes(i)));
    model.fermentationJoindown(i) = sum(ismember(fermentationJoindown,model.genes(i)));
end
%% Look at pathways
Pathway={};
for i=4:size(model,2)
    Pathway{1,i-3}=varfun(@sum,model(model{:,i}>0,:),'InputVariables',model.Properties.VariableNames{i},...
       'GroupingVariables','subSystems');
end