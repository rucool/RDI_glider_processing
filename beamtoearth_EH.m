% THis beam to earth function looks a lot more reliable than the other one
% we have. It uses all three (pitch, heading, roll) matricies to convert
% the data

function data=beamtoearth(data,cfg) %,cut) %Commented out by SC
%X	c*a	-c*a	0	0
%Y	0	0	-c*a	c*a
%Z	B	b	B	b
%E	D	d	-d	-d%
%Where:
%c =	+1 for a convex transducer head, -1 for concave
%a =	1/ [2 sin()] = 1.4619 for  =  20°, 1.0000 for  =  30°
%b =	1/ [4 cos()] = 0.2660 for  =  20°, 0.2887 for  =  30°
%d =	 = 1.0337 for  =  20°, 0.7071 for  =  30°
%© RD Instruments Version 1.0 (Jan. 11, 1999)

%data=odata.data;
%cfg=odata.cfg;
degrad=pi/180;
a=1/(2*sin(cfg.bmangle*degrad));
b=1/(4*cos(cfg.bmangle*degrad));
d=a/sqrt(2);

magdec=-12.0;%Magnetic decliniaton

xformat=zeros(4,4);
xformat(1,1)=a;
xformat(1,2)=-a;
xformat(2,3)=-a;
xformat(2,4)=a;
xformat(3,1:4)=b;
xformat(4,1:2)=d;
xformat(4,3:4)=-d;
xformat=xformat';

for i=1:length(data.time)
    
    
    %warning off
    %data=mapdepthcell(data,cfg);
    %warning on
    
    
    ubeam=[data.u1(:,i) data.u2(:,i) data.u3(:,i) data.u4(:,i)];
    
    
    
    % Transform to x,y,z coords using the matrix for the instrument
    uinstr=(ubeam*xformat); 
    
    sh=sin((data.heading(i)+magdec)*degrad);
    ch=cos((data.heading(i)+magdec)*degrad);
    sr=sin((data.roll(i)+180)*degrad); % Why +180
    cr=cos((data.roll(i)+180)*degrad);
    p=atan(tan(data.pitch(i)*degrad)*cr);
    sp=sin(p);
    cp=cos(p);
    uinstr(:,1)=uinstr(:,1);
    uinstr(:,2)=uinstr(:,2);
    uinstr(:,3)=uinstr(:,3);
    
    
    xform1=[cr 0 sr;...
        0 1 0;...
        -sr 0 cr];
    uenu=(uinstr(:,1:3)*xform1');
    
    xform2=[1 0 0;...
        0 cp sp;...
        0 -sp cp];
    
    uenu=(uenu(:,1:3)*xform2');
    
    xform3=[ch sh 0;...
        -sh ch 0;...
        0 0 1];
    
    uenu=(uenu(:,1:3)*xform3');
    
    data.VE(:,i)=uenu(:,1);
    data.VN(:,i)=uenu(:,2);
    data.VU(:,i)=uenu(:,3);
    data.VER(:,i)=uinstr(:,4);
    
    
    % Commented out by SC
%     Id=find(data.bins>cut*data.pressure(i));
%     
%     data.VE(Id,i)=nan;
%     data.VN(Id,i)=nan;
%     data.VU(Id,i)=nan;
%     data.VER(Id,i)=nan;
%     data.ei1(Id,i)=nan;
%     data.ei2(Id,i)=nan;
%     data.ei3(Id,i)=nan;
%     data.ei4(Id,i)=nan;
    
    
    
end
%odata.data=data;
end