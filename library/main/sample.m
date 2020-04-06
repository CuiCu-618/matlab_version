function [s,wt] = sample(element,s)
% !
% ! This subroutine returns the local coordinates and weighting coefficients
% ! of the integrating points.
% !
root3 = 1/sqrt(3);
r15 = 0.2/sqrt(15);
nip = size(s,1);
w = [5/9,8/9,5/9];
v = [5/9*w,8/9*w,5/9*w];
switch element
    case "line"
        switch nip
            case 1
                s(1,1) = 0;
                wt(1) = 2;
            case 2
                s(1,1)=-0.577350269189626;
                s(2,1)= 0.577350269189626;
                wt(1) = 1.000000000000000;
                wt(2) = 1.000000000000000;
            case 3
                s(1,1)=-0.774596669241484;
                s(2,1)= 0.000000000000000;
                s(3,1)= 0.774596669241484;
                wt(1) = 0.555555555555556;
                wt(2) = 0.888888888888889;
                wt(3) = 0.555555555555556;
            case 4
                s(1,1)=-0.861136311594053;
                s(2,1)=-0.339981043584856;
                s(3,1)= 0.339981043584856;
                s(4,1)= 0.861136311594053;
                wt(1) = 0.347854845137454;
                wt(2) = 0.652145154862546;
                wt(3) = 0.652145154862546;
                wt(4) = 0.347854845137454;
            case 5
                s(1,1)=-0.906179845938664;
                s(2,1)=-0.538469310105683;
                s(3,1)= 0.000000000000000;
                s(4,1)= 0.538469310105683;
                s(5,1)= 0.906179845938664;
                wt(1) = 0.236926885056189;
                wt(2) = 0.478628670499366;
                wt(3) = 0.568888888888889;
                wt(4) = 0.478628670499366;
                wt(5) = 0.236926885056189;
            case 6
                s(1,1)=-0.932469514203152;
                s(2,1)=-0.661209386466265;
                s(3,1)=-0.238619186083197;
                s(4,1)= 0.238619186083197;
                s(5,1)= 0.661209386466265;
                s(6,1)= 0.932469514203152;
                wt(1) = 0.171324492379170;
                wt(2) = 0.360761573048139;
                wt(3) = 0.467913934572691;
                wt(4) = 0.467913934572691;
                wt(5) = 0.360761573048139;
                wt(6) = 0.171324492379170;
            case 7
                s(1,1)=-0.9491079123427585245261897;
                s(2,1)=-0.7415311855993944398638648;
                s(3,1)=-0.4058451513773971669066064;
                s(4,1)= 0.000000000000000;
                s(5,1)= 0.4058451513773971669066064;
                s(6,1)= 0.7415311855993944398638648;
                s(7,1)= 0.9491079123427585245261897;
                wt(1) = 0.1294849661688696932706114;
                wt(2) = 0.2797053914892766679014678;
                wt(3) = 0.3818300505051189449503698;
                wt(4) = 0.4179591836734693877551020;
                wt(5) = 0.3818300505051189449503698;
                wt(6) = 0.2797053914892766679014678;
                wt(7) = 0.1294849661688696932706114;
            case 8
                s(1,1)=-0.9602898564975362316835609;
                s(2,1)=-0.7966664774136267395915539;
                s(3,1)=-0.5255324099163289858177390;
                s(4,1)=-0.1834346424956498049394761;
                s(5,1)= 0.1834346424956498049394761;
                s(6,1)= 0.5255324099163289858177390;
                s(7,1)= 0.7966664774136267395915539;
                s(8,1)= 0.9602898564975362316835609;
                wt(1) = 0.1012285362903762591525314;
                wt(2) = 0.2223810344533744705443560;
                wt(3) = 0.3137066458778872873379622;
                wt(4) = 0.3626837833783619829651504;
                wt(5) = 0.3626837833783619829651504;
                wt(6) = 0.3137066458778872873379622;
                wt(7) = 0.2223810344533744705443560;
                wt(8) = 0.1012285362903762591525314;
            case 9
                s(1,1)=-0.9681602395076260898355762;
                s(2,1)=-0.8360311073266357942994298;    
                s(3,1)=-0.6133714327005903973087020;
                s(4,1)=-0.3242534234038089290385380;    
                s(5,1)= 0.000000000000000;                            
                s(6,1)= 0.3242534234038089290385380;                            
                s(7,1)= 0.6133714327005903973087020;                            
                s(8,1)= 0.8360311073266357942994298;                            
                s(9,1)= 0.9681602395076260898355762;                            
                wt(1) = 0.0812743883615744119718922;                            
                wt(2) = 0.1806481606948574040584720;                            
                wt(3) = 0.2606106964029354623187429;                            
                wt(4) = 0.3123470770400028400686304;                            
                wt(5) = 0.3302393550012597631645251;                            
                wt(6) = 0.3123470770400028400686304;                            
                wt(7) = 0.2606106964029354623187429;                            
                wt(8) = 0.1806481606948574040584720;                            
                wt(9) = 0.0812743883615744119718922;                            
            case 10
                s(1,1)=-0.9739065285171717200779640;            
                s(2,1)=-0.8650633666889845107320967; 
                s(3,1)=-0.6794095682990244062343274; 
                s(4,1)=-0.4333953941292471907992659; 
                s(5,1)=-0.1488743389816312108848260; 
                s(6,1)= 0.1488743389816312108848260; 
                s(7,1)= 0.4333953941292471907992659; 
                s(8,1)= 0.6794095682990244062343274; 
                s(9,1)= 0.8650633666889845107320967; 
                s(10,1)= 0.9739065285171717200779640; 
                wt(1) = 0.0666713443086881375935688;                     
                wt(2) = 0.1494513491505805931457763;                     
                wt(3) = 0.2190863625159820439955349;                     
                wt(4) = 0.2692667193099963550912269;                     
                wt(5) = 0.2955242247147528701738930;                     
                wt(6) = 0.2955242247147528701738930;                      
                wt(7) = 0.2692667193099963550912269;                     
                wt(8) = 0.2190863625159820439955349;                     
                wt(9) = 0.1494513491505805931457763;                     
                wt(10) = 0.0666713443086881375935688;
            otherwise
                fprintf("Wrong number of integrating points for a line")
        end
    case "triangle"
        switch nip
            case 1
                s(1,1)= 0.333333333333333;
                s(1,2)= 0.333333333333333;
                wt(1) = 0.500000000000000;
            case 3
                s(1,1)= 0.500000000000000;
                s(1,2)= 0.500000000000000;
                s(2,1)= 0.500000000000000;
                s(2,2)= 0.000000000000000;
                s(3,1)= 0.000000000000000;
                s(3,2)= 0.500000000000000;
                wt(1:3)=0.333333333333333;
                wt=0.5*wt;
            case 4
                s(1,1)= 0.6;
                s(1,2)= 0.2;
                s(2,1)= 0.2;
                s(2,2)= 0.6;
                s(3,1)= 0.2;
                s(3,2)= 0.2;
                s(4,1)= 0.333333333333333;
                s(4,2)= 0.333333333333333;
                wt(1:3)= 0.520833333333333;
                wt(4)=  -0.5625;
                wt=0.5*wt;
            case 6
                s(1,1)= 0.816847572980459;
                s(1,2)= 0.091576213509771;
                s(2,1)= 0.091576213509771;
                s(2,2)= 0.816847572980459;
                s(3,1)= 0.091576213509771;
                s(3,2)= 0.091576213509771;
                s(4,1)= 0.108103018168070;
                s(4,2)= 0.445948490915965;
                s(5,1)= 0.445948490915965;
                s(5,2)= 0.108103018168070;
                s(6,1)= 0.445948490915965;
                s(6,2)= 0.445948490915965;
                wt(1:3)=0.109951743655322;
                wt(4:6)=0.223381589678011;
                wt=0.5*wt;
            case 7
                s(1,1)= 0.333333333333333;
                s(1,2)= 0.333333333333333;
                s(2,1)= 0.797426985353087;
                s(2,2)= 0.101286507323456;
                s(3,1)= 0.101286507323456;
                s(3,2)= 0.797426985353087;
                s(4,1)= 0.101286507323456;
                s(4,2)= 0.101286507323456;
                s(5,1)= 0.470142064105115;
                s(5,2)= 0.059715871789770;
                s(6,1)= 0.059715871789770;
                s(6,2)= 0.470142064105115;
                s(7,1)= 0.470142064105115;
                s(7,2)= 0.470142064105115;
                wt(1) = 0.225000000000000;
                wt(2:4)=0.125939180544827;
                wt(5:7)=0.132394152788506;
                wt=0.5*wt;
            case 12
                s(1,1)= 0.873821971016996;
                s(1,2)= 0.063089014491502;
                s(2,1)= 0.063089014491502;
                s(2,2)= 0.873821971016996;
                s(3,1)= 0.063089014491502;
                s(3,2)= 0.063089014491502;
                s(4,1)= 0.501426509658179;
                s(4,2)= 0.249286745170910;
                s(5,1)= 0.249286745170910;
                s(5,2)= 0.501426509658179;
                s(6,1)= 0.249286745170910;
                s(6,2)= 0.249286745170910;
                s(7,1) =0.053145049844817;
                s(7,2) =0.310352451033784;
                s(8,1) =0.310352451033784;
                s(8,2) =0.053145049844817;
                s(9,1) =0.053145049844817;
                s(9,2) =0.636502499121398;
                s(10,1)=0.310352451033784;
                s(10,2)=0.636502499121398;
                s(11,1)=0.636502499121398;
                s(11,2)=0.053145049844817;
                s(12,1)=0.636502499121398;
                s(12,2)=0.310352451033784;
                wt(1:3)=0.050844906370207;
                wt(4:6)=0.116786275726379;
                wt(7:12)=0.082851075618374;
                wt=0.5*wt;
            case 16
                s(1,1)=0.333333333333333;
                s(1,2)=0.333333333333333;
                s(2,1)=0.658861384496478;
                s(2,2)=0.170569307751761;
                s(3,1)=0.170569307751761;
                s(3,2)=0.658861384496478;
                s(4,1)=0.170569307751761;
                s(4,2)=0.170569307751761;
                s(5,1)=0.898905543365938;
                s(5,2)=0.050547228317031;
                s(6,1)=0.050547228317031;
                s(6,2)=0.898905543365938;
                s(7,1)=0.050547228317031;
                s(7,2)=0.050547228317031;
                s(8,1)=0.081414823414554;
                s(8,2)=0.459292588292723;
                s(9,1)=0.459292588292723;
                s(9,2)=0.081414823414554;
                s(10,1)=0.459292588292723;
                s(10,2)=0.459292588292723;
                s(11,1)=0.008394777409958;
                s(11,2)=0.263112829634638;
                s(12,1)=0.008394777409958;
                s(12,2)=0.728492392955404;
                s(13,1)=0.263112829634638;
                s(13,2)=0.008394777409958;
                s(14,1)=0.263112829634638;
                s(14,2)=0.728492392955404;
                s(15,1)=0.728492392955404;
                s(15,2)=0.008394777409958;
                s(16,1)=0.728492392955404;
                s(16,2)=0.263112829634638;
                wt(1)=0.144315607677787;
                wt(2:4)=0.103217370534718;
                wt(5:7)=0.032458497623198;
                wt(8:10)=0.095091634267284;
                wt(11:16)=0.027230314174435;
                wt=0.5*wt;
            otherwise
                fprintf("Wrong number of integrating points for a triangle")                
        end
    case "quadrilateral"
        switch nip
            case 1
                s(1,1)=0.0;
                s(1,2)=0.0;
                wt(1)=4.0;
            case 4
                s(1,1)=-root3;
                s(1,2)= root3;
                s(2,1)= root3;
                s(2,2)= root3;
                s(3,1)=-root3;
                s(3,2)=-root3;
                s(4,1)= root3;
                s(4,2)=-root3;
                wt=1.0;
            case 9
                s(1:7:3,1)=-r15;
                s(2:8:3,1)=0.0;
                s(3:9:3,1)=r15;
                s(1:3,2)  =r15;
                s(4:6,2)  =0.0;
                s(7:9,2)  =-r15;
                wt= v;
            case 16
                s(1:13:4,1)=-0.861136311594053;
                s(2:14:4,1)=-0.339981043584856;
                s(3:15:4,1)= 0.339981043584856;
                s(4:16:4,1)= 0.861136311594053;
                s(1:4,2)   = 0.861136311594053;
                s(5:8,2)   = 0.339981043584856;
                s(9:12,2)  =-0.339981043584856;
                s(13:16,2) =-0.861136311594053;
                wt(1)      = 0.121002993285602;
                wt(4)      = wt(1);
                wt(13)     = wt(1);
                wt(16)     = wt(1);
                wt(2)      = 0.226851851851852;
                wt(3)      = wt(2);
                wt(5)      = wt(2);
                wt(8)      = wt(2);
                wt(9)      = wt(2);
                wt(12)     = wt(2);
                wt(14)     = wt(2);
                wt(15)     = wt(2);
                wt(6)      = 0.425293303010694;
                wt(7)      = wt(6);
                wt(10)     = wt(6);
                wt(11)     = wt(6);
            case 25
                s(1:21:5,1)= 0.906179845938664;
                s(2:22:5,1)= 0.538469310105683;
                s(3:23:5,1)= 0.0;
                s(4:24:5,1)=-0.538469310105683;
                s(5:25:5,1)=-0.906179845938664;
                s( 1: 5,2) = 0.906179845938664;
                s( 6:10,2) = 0.538469310105683;
                s(11:15,2) = 0.0;
                s(16:20,2) =-0.538469310105683;
                s(21:25,2) =-0.906179845938664;
                wt(1) =0.056134348862429;
                wt(2) =0.113400000000000;
                wt(3) =0.134785072387521;
                wt(4) =0.113400000000000;
                wt(5) =0.056134348862429;
                wt(6) =0.113400000000000;
                wt(7) =0.229085404223991;
                wt(8) =0.272286532550750;
                wt(9) =0.229085404223991;
                wt(10)=0.113400000000000;
                wt(11)=0.134785072387521;
                wt(12)=0.272286532550750;
                wt(13)=0.323634567901235;
                wt(14)=0.272286532550750;
                wt(15)=0.134785072387521;
                wt(16)=0.113400000000000;
                wt(17)=0.229085404223991;
                wt(18)=0.272286532550750;
                wt(19)=0.229085404223991;
                wt(20)=0.113400000000000;
                wt(21)=0.056134348862429;
                wt(22)=0.113400000000000;
                wt(23)=0.134785072387521;
                wt(24)=0.113400000000000;
                wt(25)=0.056134348862429;
            otherwise
                fprintf("Wrong number of integrating points for a quadrilaterals")                                
        end
    case "tetrahedron"
        % for tetrahedra weights multiplied by 1/6
        switch nip
            case 1
                s(1,1)=0.25;
                s(1,2)=0.25;
                s(1,3)=0.25;
                wt(1)=1.0/6.0;
            case 4
                s(1,1)=0.58541020;
                s(1,2)=0.13819660;
                s(1,3)=s(1,2);
                s(2,2)=s(1,1);
                s(2,3)=s(1,2);
                s(2,1)=s(1,2);
                s(3,3)=s(1,1);
                s(3,1)=s(1,2);
                s(3,2)=s(1,2);
                s(4,1)=s(1,2);
                s(4,2)=s(1,2);
                s(4,3)=s(1,2);
                wt(1:4)=0.25/6.0;
            case 5
                s(1,1)=0.25;
                s(1,2)=0.25;
                s(1,3)=0.25;
                s(2,1)=0.5;
                s(2,2)=1.0/6.0;
                s(2,3)=s(2,2);
                s(3,2)=0.5;
                s(3,3)=1.0/6.0;
                s(3,1)=s(3,3);
                s(4,3)=0.5;
                s(4,1)=1.0/6.0;
                s(4,2)=s(4,1);
                s(5,1)=1.0/6.0;
                s(5,2)=s(5,1);
                s(5,3)=s(5,1);
                wt(1)=-0.8;
                wt(2)=9.0/20.0;
                wt(3:5)=wt(2);
                wt=wt/6.0;
            otherwise
                fprintf("Wrong number of integrating points for a tetrahedron")                                                
        end
    case "hexahedron"
        switch nip
            case 1
                s(1,1:3)=0.0;
                wt(1)=8.0;
            case 8
                s(1,1)= root3;
                s(1,2)= root3;
                s(1,3)= root3;
                s(2,1)= root3;
                s(2,2)= root3;
                s(2,3)=-root3;
                s(3,1)= root3;
                s(3,2)=-root3;
                s(3,3)= root3;
                s(4,1)= root3;
                s(4,2)=-root3;
                s(4,3)=-root3;
                s(5,1)=-root3;
                s(5,2)= root3;
                s(5,3)= root3;
                s(6,1)=-root3;
                s(6,2)=-root3;
                s(6,3)= root3;
                s(7,1)=-root3;
                s(7,2)= root3;
                s(7,3)=-root3;
                s(8,1)=-root3;
                s(8,2)=-root3;
                s(8,3)=-root3;
                wt=1.0;
            case 14
                b=0.795822426;
                c=0.758786911;
                wt(1:6)=0.886426593;
                wt(7:14)=0.335180055;
                s(1,1)=-b;
                s(2,1)=b;
                s(3,2)=-b;
                s(4,2)=b;
                s(5,3)=-b;
                s(6,3)=b;
                s(7:end,:)=c;
                s(7,1)=-c;
                s(7,2)=-c;
                s(7,3)=-c;
                s(8,2)=-c;
                s(8,3)=-c;
                s(9,1)=-c;
                s(9,3)=-c;
                s(10,3)=-c;
                s(11,1)=-c;
                s(11,2)=-c;
                s(12,2)=-c;
                s(13,1)=-c;
            case 15
                b=1.0;
                c       =0.674199862;
                wt(1)   =1.564444444;
                wt(2:7) =0.355555556;
                wt(8:15)=0.537777778;
                s(2,1)=-b;
                s(3,1)=b;
                s(4,2)=-b;
                s(5,2)=b;
                s(6,3)=-b;
                s(7,3)=b;
                s(8:end,:)=c;
                s(8,1)=-c;
                s(8,2)=-c;
                s(8,3)=-c;
                s(9,2)=-c;
                s(9,3)=-c;
                s(10,1)=-c;
                s(10,3)=-c;
                s(11,3)=-c;
                s(12,1)=-c;
                s(12,2)=-c;
                s(13,2)=-c;
                s(14,1)=-c;
            case 27
                wt=[5.0/9.0*v,8.0/9.0*v,5.0/9.0*v];
                s(1:7:3,1)=-r15;
                s(2:8:3,1)=0.0;
                s(3:9:3,1)=r15;
                s(1:3,3)=r15;
                s(4:6,3)=0.0;
                s(7:9,3)=-r15;
                s(1:9,2)=-r15;
                s(10:16:3,1)=-r15;
                s(11:17:3,1)=0.0;
                s(12:18:3,1)=r15;
                s(10:12,3)=r15;
                s(13:15,3)=0.0;
                s(16:18,3)=-r15;
                s(10:18,2)=0.0;
                s(19:25:3,1)=-r15;
                s(20:26:3,1)=0.0;
                s(21:27:3,1)=r15;
                s(19:21,3)= r15;
                s(22:24,3)=0.0;
                s(25:27,3)=-r15;
                s(19:27,2)= r15;
            otherwise
                fprintf("Wrong number of integrating points for a hexahedron")
        end
    otherwise
                fprintf("not a valid element type")        
end
end