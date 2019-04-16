function [input_index, output_index, eta0, spot, rock_type, sanoffset, phioffset] = input_record(input)

%defaults
input_index = 1;
output_index = 1;
eta0 = 0.8337;
spot = [80 76];
rock_type = 'san';
sanoffset =0;
phioffset =0;


    switch input
        
  
        
        case 'Hazuki'
            input_index = 1; %data location (number)
            output_index = 1; %where to put result
            eta0 = 0.28;  %rocking width for spot in degrees
            spot = [80 80];  %Spot coordinates
            sanoffset = -0.2;  %san misalignment
            phioffset = 0.0;  %phi misalignment
            rock_type = 'san';
        

        
        case 'Heidi omega 175mT'
            input_index = 1; %data location (number)
            output_index = 1; %where to put result
            eta0 = 0.09;  %rocking width for spot in degrees
            spot = [91 59];  %Spot coordinates
            sanoffset = -89.9;  %san misalignment
            phioffset =0.04;  %phi misalignment
            rock_type = 'san';
        
        case 'PdBi2 10 deg. off'
            
            eta0 = 0.113   %rocking fwhm in deg
            spot = [88 60]  %[x0 y0]
            sanoffset = -90.11;
            phioffset = -0.13;
            rock_type = 'phi';
            
        case 'Ru7B3_0.2T'
            input_index = 1;
            output_index = 1;
            eta0 =  0.3;
            spot = [178 67];
            rock_type = 'san';
            sanoffset =0;
            phioffset =0;
        
        case 'CeCu2Si2_0.8T' 
            input_index = 4;
            output_index = 4;
            eta0 =  0.1;
            spot = [80 90];
            rock_type = 'phi';
            sanoffset =0.0083498;
            phioffset =0.0048556;
        case 'CeCu2Si2_2T' 
            input_index = 3;
            output_index = 3;
            eta0 =  0.111856;
            spot = [60 103];
            rock_type = 'phi';
            sanoffset =0.0083498;
            phioffset =0.0048556;

        case 'CeCu2Si2_1.8T' 
            input_index = 2;
            output_index = 2;
            eta0 =  0.18535;
            spot = [63 100];
            rock_type = 'phi';
            sanoffset =0.0083498;
            phioffset =0.0048556;


        case 'CeCu2Si2_1.5T' 
            input_index = 1;
            output_index = 1;
            eta0 =  0.1341;
            spot = [67 97];
            rock_type = 'phi';
            sanoffset =0.0083498;
            phioffset =0.0048556;

        case 'CaYBCO 7T' 
            input_index = 2;
            output_index = 2;
            eta0 =  1.1;
            spot = [180 38];
            rock_type = 'phi';
            sanoffset =-0.0044212;
            phioffset =0.017692;

        case 'CaYBCO 16.4T' 
            input_index = 1;
            output_index = 1;
            eta0 =  0.4464;
            spot = [214 106];
            rock_type = 'phi';
            sanoffset =-0.013159;
            phioffset =-0.002503;
      
        case 'BFAP 5T san'
            input_index = 1;
            output_index = 1;
            %eta0 = 0.436833;
            eta0 = 0.440735;
            %eta0=0.6
            spot = [91 65];
            rock_type = 'san';
            sanoffset =0.0260665;
            phioffset =-0.0371592;
            %sanoffset =0.0258007;
            %phioffset =-0.0384957;
            
       case 'BFAP PSI 1T sanphi'
            input_index = 2;
            output_index = 2;
            eta0 = 1;
            spot = [80 76];
            rock_type = 'san';
            sanoffset =0.0421451;
            phioffset =0.0315147;
            
         case 'BFAP 0.7T sanphi'
            input_index = 1;
            output_index = 1;
            eta0 = 1.16;
            spot = [95 91];
            rock_type = 'san';
            sanoffset =0.0421451;
            phioffset =0.0315147;
 
        case 'BFAP 10T san'
            input_index = 1;
            output_index = 1;
            eta0 = 0.480425;
            spot = [102 64];
            rock_type = 'san';
            sanoffset =0.0258007;
            phioffset =-0.0384957;
         case 'BFAP 7T san'
            input_index = 1;
            output_index = 1;
            eta0 = 0.473978;
            spot = [96 64];
            rock_type = 'san';
            sanoffset =0.0258007;
            phioffset =-0.0384957;

        case 'BFAP 3T san'
            input_index = 1;
            output_index = 1;
            eta0 = 0.833231;
            spot = [84 97];
            rock_type = 'san';
            sanoffset =0.0258007;
            phioffset =-0.0384957;

        case 'BFAP 1T san'
            input_index = 1;
            output_index = 1;
            eta0 = 1.30019;
            spot = [52 81];
            rock_type = 'san';
            sanoffset =0.0258007;
            phioffset =-0.0384957;
         case '124 0.4T'
              %Sector boxes 13 22 90 40
              %             20 30 160 25
              %             20 30 200 25
              %             13 22 270 40
              %             20 30 -20 25
              %             20 30 20 25
            input_index = 1;
            output_index = 1;
            eta0 =  2.7400;
            spot = [72 85];
            rock_type = 'phi';
            sanoffset = 0;
            phioffset =  0.11052;
           case '124 0.8T'
              %Sector boxes 30 40 90 30
              %             40 55 158 25
              %             40 55 204 25
              %             30 40 270 30
              %             40 55 -24 25
              %             40 55 24 25
            input_index = 2;
            output_index = 2;
            eta0 =  1.4865;
            spot = [49 105];
            rock_type = 'phi';
            sanoffset = -0.0517081;
            phioffset =  0.13668;
          case '124 4T'
            input_index = 3;
            output_index = 3;
            eta0 =  1.0607;
            spot = [78 91];
            rock_type = 'phi';
            sanoffset = -0.012863;
            phioffset =  0.17143;
        case '124 6T'
            input_index = 4;
            output_index = 4;
            eta0 =  0.707736;
            spot = [81 96];
            rock_type = 'phi';
            sanoffset = 0.0529769;
            phioffset = 0.0901328;
        case 'FRM2 Nb'
            input_index = 1;
            output_index = 1;
            eta0 =  179%0.147501;
            spot = [44 62];
            rock_type = 'san';
            sanoffset =0%-0.270658;
            phioffset =0%0.18446;
        case '124 1a'
            input_index = 1;
            output_index = 1;
            eta0 =  0.974545;
            spot = [82 96];
            rock_type = 'phi';
            sanoffset = -0.076443;
            phioffset = 0.149058;
        case 'CeCoIn5 HZB'
            input_index = 3;
            output_index = 3;
            eta0 = 0.3;
            spot = [49 47];
            rock_type = 'phi';
            sanoffset =-1.3845;
            phioffset = -0.36352;
          
        case 'BiPd 100 600Gphi'
            input_index = 3;
            output_index = 3;
            eta0 = 1.72005;
            spot = [69 68];
            rock_type = 'phi';
            sanoffset =-0.0725148;
            phioffset =-0.0350695;

        case 'BiPd 100 250Gphi'
            input_index = 2;
            output_index = 2;
            eta0 = 1.5788;
            spot = [69 68];
            rock_type = 'phi';
            sanoffset =-0.034727;
            phioffset =0.0705434;
 
        case 'BiPd 450Gphi'
            input_index = 1;
            output_index = 1;
            eta0 = 0.294763;
            spot = [40 64];
            rock_type = 'phi';
            sanoffset =-0.034727;
            phioffset =0.0705434;
        case 'BiPd 400Gphi'
            input_index = 1;
            output_index = 1;
            eta0 = 0.248735;
            spot = [40 64];
            rock_type = 'phi';
            sanoffset =-0.0324092;
            phioffset =0.0599814;
        case 'BiPd 200Gphi'
            input_index = 3;
            output_index = 3;
            eta0 = 1.3;
            spot = [40 64];
            rock_type = 'phi';
            sanoffset =0;
            phioffset =0;
       case 'BiPd 200Gphi b'
            input_index = 2;
            output_index = 2;
            eta0 = 1.3;
            spot = [63 72];
            rock_type = 'phi';
            sanoffset =0;
            phioffset =0;
         case 'BiPd 250Gsan a'
            input_index = 1;
            output_index = 1;
            eta0 = 1.06;
            spot = [159 65];
            rock_type = 'phi';
            sanoffset =0;
            phioffset =0;
        case 'BiPd 250Gphi a'
            input_index = 1;
            output_index = 3;
            eta0 = 1.06;
            spot = [155 69];
            rock_type = 'phi';
            sanoffset =0;
            phioffset =0;
        
        case 'BFAP 0.7T san'
            input_index = 1;
            output_index = 1;
            eta0 = 2.3098;
            spot = [134 175];
            rock_type = 'san';
            sanoffset =0.565643;
            phioffset =0;
            
            
        case 'BFAP 0.7T phi'
            fit=0;
            input_index = 2;
            output_index = 4;
            eta0 = 2.3098;
            spot = [134 175];
            rock_type = 'san';
            sanoffset =0.565643;
            phioffset =0;
           
        case 'BFAP 12T'
            eta0 = 1.2221
            spot = [20 64]
            rock_type = 'san'
   
        case 'Ca doped YBCO 16.4T phi'
            input_index = 2; %data location
            output_index = 2;
            eta0 = 0.79139;
            spot = [181 105];
            sanoffset =0;
            phioffset =0;

        case 'BFAP 0.1T san'
            %FG HiResSANS_exp12_scan0040_
            %BG HiResSANS_exp12_scan0042_
            input_index = 1; %data location
            output_index = 1;

            spot = [115 148];
            eta0 =   1;
            rock_type = 'phi'
            sanoffset = 0.552682;
            phioffset = 0.011199;

            
        case 'BiPd'
            eta0 = 1.0201
            %eta0 = 0.1779
            %eta0 = 0.82194
            % fitted eta0 = 0.6934
            spot = [103 74]
            rock_type = 'san'
        case 'Nbphi'
            input_index = 1; %data location
            output_index = 1;
            %eta0 = 0.3782;
            eta0 =   0.55663;
            %eta0 = 1.7664
            spot = [115 148];
            eta0 =   0.55663;
            rock_type = 'phi'
            sanoffset = 0.552682;
            phioffset = 0.011199;
          case 'Nbsan'
            input_index = 3; %data location
            output_index = 3;
            spot = [116 148];
            eta0 =  0.699838;
            rock_type = 'san'
            sanoffset = 0.565703;
            phioffset = -0.0180286;
           
            
        case 'KFA_1'
            eta0 = 0.75;
            spot = [45 91]
            
        case 'CeCoIn5 4.55T'
            eta0 = 0.1891;
            spot = [88 74]
            
        case 'CeCoIn5 9.5T'
            eta0 = 0.32263;
            spot = [82 87]
            
        case 'BFAP 12T'
            %eta0 =  0.57688;
            eta0 = 0.57688;
            spot = [110 65]
        case 'YBCO 10T 30deg'
            eta0 = 0.82
            %eta0 = 0.6136
            spot = [90 88]
        case 'YBCO6.5'
            eta0 = 1.038;
            spot = [84 66]
            
        case 'YBCO 30deg'
            eta0 = 1.038;
            spot = [79 80]
        case 'LSCO 1T'
            eta0 = 1.5247;
            spot=  [58 105];
        case '3T'
            eta0 = 1.5247;
            spot=  [84 97];
        case '5T'
            % 5T:
            %eta0 =  0.4991
            eta0 =   0.4232 %fitted
            %eta0 = 0.7998
            %eta0 = 0.63621
            %eta0 = 0.3494
            %eta0 =  0.4211;  %rocking curve width degrees for on-axis spot get from fit
            sanoffset = 0.026058
            phioffset = -0.037127
            spot=  [37 64]; %spot position of on-axis spot in pixels
            rock_type = 'san'
        case '10T'
            %10T 
            eta0 = 0.4844
            %eta0 = 0.68752 %rocking curve FWHM degrees for on-axis spot get from fit
            spot = [103 64]  %spot position
            %     case '5T'
            %         % 5T:
            %         eta0 =  0.4211  %rocking curve width degrees for on-axis spot get from fit
            %         x0 = 91.4765  %spot position of on-axis spot in pixels
        case '1T'
            %1 T
            eta0 =  1.0936;  %rocking curve width degrees for on-axis spot get from fit
            x0 = 78  %spot position of on-axis spot in pixels
        case '7T'
            %7T
            %eta0 =  0.48589
            eta0 = 0.6670
            %eta0 = 6.8949
            spot = [96 64]
            rock_type = 'san'
        case '30degsan'
            %0.2T 30 deg
            %san rock EB 582 FG 583:594,596:607 BG 633>641,643>650
            dataset = 1
            eta0 =   0.4012;
            spot = [39 52]
            sanoffset = 0;
            phioffset = 0;
            input_index = 1
            output_index = 1
            rock_type = 'san'
        case '30degphi'
            %phi rock EB 582 FG 608:617,619:632 BG 658>666,668>678, alpha = 90
            fixed = [0 1 1];
            eta0 =   0.382481;
            spot = [63 86];
            rock_type = 'phi'
            input_index = 2
            output_index = 3
        case '45degsan'
            dataset =1
            eta0= 0.36263
            x0= 85
            input_index = 1
            output_index = 3
        case '45degphi'
            dataset = 2
            eta0= 0.36263
            x0= 85
            alpha = 90;
            input_index = 2
            output_index = 5
            ignorant = 0;
            priorset = 1;
        case '75degsan'
            dataset =1
            %eta0= 0.36263
            
            x0= 85
            input_index = 1
        case '75degphi'
            dataset = 2
            %eta0= 0.36263
            x0= 85
            alpha = 90;
            input_index = 2
        case '15degsan'
            input_index=1
            dataset =1
            eta0 = 0.48787
            x0=73
        case '15degphi'
            input_index=2
            dataset =2
            eta0 = 0.48787
            x0=73
            alpha=90
            
    end


end