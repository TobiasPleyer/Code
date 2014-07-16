
in = dlmread('../Daten/26_AIR_FROG_31.0W_RTT=7.415us_Ip=12.2A.bin.Speck.dat');
c = 299792458;
oo = 2*pi*c ./ in(:,1).* 1e9;
oi = in(:,2);
op = in(:,3);

in = dlmread('filtered_original.txt');
fo = in(:,1);
fi = in(:,2);
fp = in(:,3);

in = dlmread('filtered_super_smoothed.txt');
so = in(:,1);
si = in(:,2);
sp = in(:,3);

figure(1)
plot(oo,oi,'LineWidth',3)
hold on
plot(fo,fi,'color','red','LineWidth',2)
plot(so,si,'color','green')
hold off

figure(2)
plot(oo,op,'LineWidth',3)
hold on
plot(fo,fp,'color','red','LineWidth',2)
plot(so,sp,'color','green')
hold off