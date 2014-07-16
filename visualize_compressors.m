GDD = zeros(2,50);
TOD = zeros(2,50);

grating_step_size   = 0.15/50;
prism_step_size     = 4.5/50;
grating_start_point = 0.05; % i.e we go from 5cm to 20cm
prism_start_point   = 0.5;  % i.e we go from 50cm to 5m

grating_range = grating_start_point:grating_step_size:0.2;
prism_range   = prism_start_point:prism_step_size:5;

for i=0:49
    [GDD(1,i),TOD(1,i)] = Grating(grating_start_point + i*grating_step_size);
    [GDD(2,i),TOD(2,i)] = PrismS('SF11',prism_start_point + i*prism_step_size,0.05);  
end

subplot(2,2,1)
plot(grating_range,GDD(1,:))
subplot(2,2,2)
plot(grating_range,TOD(1,:))
subplot(2,2,3)
plot(prism_range,GDD(1,:))
subplot(2,2,4)
plot(prism_range,TOD(1,:))