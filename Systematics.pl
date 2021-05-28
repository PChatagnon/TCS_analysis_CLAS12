#!/apps/bin/perl



$PI = 3.14159265359;
#system( "root -l \"Plot.C(\\\"BSAversusM\\\",\\\"BSAversusM\\\",\\\"BSAversusM\\\",4,\\\"BSA\\\",\\\"M (GeV)\\\",\\\"BSA\\\",\\\"M (GeV)\\\",-0.3,1,\\\"<E#gamma>=7.3 GeV,<-t>=0.37 GeV^{2}\\\")\"");
system( "root -l \"Plot.C(\\\"BSAversust\\\",\\\"BSAversust\\\",\\\"BSAversust\\\",4,\\\"BSA\\\",\\\"-t (GeV^{2})\\\",\\\"BSA\\\",\\\"-t (GeV^{2})\\\",-0.3,1.,\\\"<E#gamma>=7.3 GeV,<M>=1.8 GeV\\\")\"");
#system( "root -l \"Plot.C(\\\"BSAversusXi\\\",\\\"BSAversusXi\\\",\\\"BSAversusXi\\\",3,\\\"BSA\\\",\\\"#xi\\\",\\\"BSA\\\",\\\"xi\\\",-0.3,1,\\\"<E#gamma>=7.3 GeV,<M>=1.8 GeV,<-t>=0.37 GeV^{2}\\\")\"");
#system( "root -l \"Plot.C(\\\"AFBXi\\\",\\\"AFBXi\\\",\\\"AFBXi\\\",3,\\\"A_{FB}\\\",\\\"#xi\\\",\\\"A_{FB}\\\",\\\"xi\\\",-0.3,1,\\\"<E#gamma>=7.2 GeV,<M>=1.81 GeV,<-t>=0.0.36 GeV^{2}\\\")\"");
system( "root -l \"Plot.C(\\\"AFBt2\\\",\\\"AFBt2\\\",\\\"AFBt2\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.3,1.,\\\"<E#gamma>=8.13 GeV,<M>=2.25 GeV\\\")\"");
#system( "root -l \"Plot.C(\\\"AFBt1\\\",\\\"AFBt1\\\",\\\"AFBt1\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.3,1,\\\"<E#gamma>=7.03 GeV,<M>=1.71 GeV\\\")\"");
system( "root -l \"Plot.C(\\\"AFBt\\\",\\\"AFBt\\\",\\\"AFBt\\\",4,\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",\\\"A_{FB}\\\",\\\"-t (GeV^{2})\\\",-0.3,1,\\\"<E#gamma>=7.2 GeV,<M>=1.81 GeV\\\")\"");
#system( "root -l \"Plot.C(\\\"AFBQ2\\\",\\\"AFBQ2\\\",\\\"AFBQ2\\\",4,\\\"A_{FB}\\\",\\\"M (GeV)\\\",\\\"A_{FB}\\\",\\\"M (GeV)\\\",-0.5,1.,\\\"<E#gamma>=7.2 GeV,<-t>=0.36 GeV^{2}\\\")\"");
#system( "root -l \"Plot.C(\\\"AFBEg\\\",\\\"AFBEg\\\",\\\"AFBEg\\\",3,\\\"A_{FB}\\\",\\\"E#gamma (GeV)\\\",\\\"A_{FB}\\\",\\\"E#gamma (GeV)\\\",-0.3,1,\\\"<-t>=0.36 GeV^{2},<M>=1.81 GeV\\\")\"");
#system( "root -l \"Plot.C(\\\"RratioXi\\\",\\\"RratioXi\\\",\\\"RratioXi\\\",3,\\\"R'\\\",\\\"#xi\\\",\\\"R'\\\",\\\"xi\\\",-0.3,0.5,\\\"<E#gamma>=7.2 GeV,<M>=1.8 GeV, <-t>=0.37 GeV^{2}\\\")\"");
#system( "root -l \"Plot.C(\\\"RRatiovst\\\",\\\"RRatiovst\\\",\\\"RRatiovst\\\",4,\\\"R'\\\",\\\"-t (GeV^{2})\\\",\\\"R'\\\",\\\"-t (GeV^{2})\\\",-0.3,0.5,\\\"<E#gamma>=7.3 GeV,<M>=1.8 GeV\\\")\"");
#system( "root -l \"Plot.C(\\\"RRatiovstHighMass\\\",\\\"RRatiovstHighMass\\\",\\\"RRatiovstHighMass\\\",4,\\\"R'\\\",\\\"-t (GeV^{2})\\\",\\\"R'\\\",\\\"-t (GeV^{2})\\\",-0.3,0.5,\\\"<E#gamma>=8.2 GeV,<M>=2.23 GeV\\\")\"");
#$cmdVGG1 = "root -l FromVGG.C out1.txt out2.txt Assym${t2}.root";
#		system($cmdVGG1);



