#!/apps/bin/perl



$PI = 3.14159265359;
system( "root -l \"PlotNoSyst.C(\\\"AFBt2\\\",\\\"AFBt2\\\",\\\"AFBt2\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.4,0.7,\\\"<E#gamma>=7.88 GeV,<M>=2.21 GeV\\\",\\\"2.21\\\")\"");
#system( "root -l \"PlotNoSyst.C(\\\"AFBt1\\\",\\\"AFBt1\\\",\\\"AFBt1\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.4,0.7,\\\"<E#gamma>=6.5 GeV,<M>=1.69 GeV\\\",\\\"1.69\\\")\"");
system( "root -l \"PlotNoSyst.C(\\\"AFBt\\\",\\\"AFBt\\\",\\\"AFBt\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.4,0.7,\\\"<E#gamma>=7.2 GeV,<M>=1.81 GeV\\\",\\\"1.79\\\")\"");

#system( "root -l \"PlotNoSyst.C(\\\"AFBQ2\\\",\\\"AFBQ2\\\",\\\"AFBQ2\\\",4,\\\"A_{FB}\\\",\\\"M (GeV)\\\",\\\"A_{FB}\\\",\\\"M (GeV)\\\",-0.55,0.6,\\\"<E#gamma>=7.2 GeV,<-t>=0.0.36 GeV^{2}\\\",\\\"AFB_M\\\")\"");

#system( "root -l \"PlotNoSyst.C(\\\"AFBEg\\\",\\\"AFBEg\\\",\\\"AFBEg\\\",3,\\\"A_{FB}\\\",\\\"E_{#gamma} (GeV)\\\",\\\"A_{FB}\\\",\\\"E_{#gamma} (GeV)\\\",-0.3,0.6,\\\"<-t>=0.36 GeV^{2},<M>=1.81 GeV\\\",\\\"AFB_Eg\\\")\"");

system( "root -l \"PlotNoSyst.C(\\\"BSAversust\\\",\\\"BSAversust\\\",\\\"BSAversust\\\",4,\\\"BSA\\\",\\\"-t (GeV^{2})\\\",\\\"BSA\\\",\\\"-t (GeV^{2})\\\",-0.3,0.6,\\\"<E#gamma>=6.8 GeV,<M>=1.79 GeV\\\",\\\"BSA_t\\\")\"");

#system( "root -l \"PlotNoSyst.C(\\\"BSAversusM\\\",\\\"BSAversusM\\\",\\\"BSAversusM\\\",4,\\\"BSA\\\",\\\"M (GeV)\\\",\\\"BSA\\\",\\\"M (GeV)\\\",-0.3,0.6,\\\"<E#gamma>=6.8 GeV,<-t>=0.33 GeV^{2}\\\",\\\"BSA_M\\\")\"");
