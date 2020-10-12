function data=mapdepthcell_EH(data,cfg)
%data=odata.data;
%cfg=odata.cfg;

% Convert vars from cfg to doubles SC
bins=((([1:double(cfg.nbins)]-1)*double(cfg.binsize))+double(cfg.bin1))/100;
data.bins=double(bins);
y=double(bins);
xp=(y*tan(20*pi/180));
xm=-(y*tan(20*pi/180));
degrad=pi/180;

for i=1:length(data.time)
    p=(data.pitch(i))*degrad;
    zp=xp+1i*y;
    zm=xm+1i*y;
    %map beams 3 and 4
    p=(data.pitch(i))*degrad;
    zp=xp+1i*y;
    zm=xm+1i*y;
    u3bins=(zp*exp(1i*p));
    u4bins=(zm*exp(1i*p));
    r=(data.roll(i))*degrad;
    u1bins=(zm*exp(1i*r));
    u2bins=(zp*exp(1i*r));
    
    %ubeam(:,3)=interp1(imag(u3bins),ubeam(:,3),bins,'linear','extrap');
    %ubeam(:,4)=interp1(imag(u4bins),ubeam(:,4),bins,'linear','extrap');
    %ubeam(:,1)=interp1(imag(u1bins),ubeam(:,1),bins,'linear','extrap');
    %ubeam(:,2)=interp1(imag(u2bins),ubeam(:,2),bins,'linear','extrap');
    warning off
    data.u3(:,i)=interp1(imag(u3bins),data.u3(:,i),bins,'linear','extrap');
    data.u4(:,i)=interp1(imag(u4bins),data.u4(:,i),bins,'linear','extrap');
    data.u1(:,i)=interp1(imag(u1bins),data.u1(:,i),bins,'linear','extrap');
    data.u2(:,i)=interp1(imag(u2bins),data.u2(:,i),bins,'linear','extrap');
    
    data.ei3(:,i)=interp1(imag(u3bins),data.ei3(:,i),bins,'linear','extrap');
    data.ei4(:,i)=interp1(imag(u4bins),data.ei4(:,i),bins,'linear','extrap');
    data.ei1(:,i)=interp1(imag(u1bins),data.ei1(:,i),bins,'linear','extrap');
    data.ei2(:,i)=interp1(imag(u2bins),data.ei2(:,i),bins,'linear','extrap');
    
    warning on
    
    
    
end

end
