[GDDg1,TODg1] = Grating(0.05);
[GDDg2,TODg2] = Grating(0.20);
[GDDp1,TODp1] = PrismS('SF11',0.50,0.05);
[GDDp2,TODp2] = PrismS('SF11',5.00,0.05);

fprintf('Bei Variation der Gitter von 5cm bis 20cm\nerhalten wir eine Range von:\nGDD: [%2.1e,%2.1e] fs^2\nTOD: [%2.1e,%2.1e] fs^3\n',GDDg1,GDDg2,TODg1,TODg2)
fprintf('\n')
fprintf('Bei Variation der Prismen von 50cm bis 5m\nund 5cm opt. Weg erhalten wir eine Range von:\nGDD: [%2.1e,%2.1e] fs^2\nTOD: [%2.1e,%2.1e] fs^3\n',GDDp1,GDDp2,TODp1,TODp2)
% GDD = zeros(2,50);
% TOD = zeros(2,50);
% 
% N = 50; % number of probe points
% 
% grating_start_point = 0.05; % i.e we start from 5cm
% prism_start_point   = 0.5;  % i.e we start from 50cm
% 
% grating_end_point = 0.2; % i.e we end at 20cm
% prism_end_point   = 5;   % i.e we end at 5m
% 
% grating_step_size = (grating_end_point-grating_start_point)/N;
% prism_step_size   = (prism_end_point-prism_start_point)/N;
% 
% grating_range = grating_start_point:grating_step_size:(grating_end_point-grating_step_size);
% prism_range   = prism_start_point:prism_step_size:(prism_end_point-prism_step_size);
% 
% for i=0:49
%     [GDD(1,i+1),TOD(1,i+1)] = Grating(grating_start_point + i*grating_step_size);
%     [GDD(2,i+1),TOD(2,i+1)] = PrismS('SF11',prism_start_point + i*prism_step_size,0.05);  
% end
% 
% subplot(2,2,1)
% plot(grating_range,GDD(1,:))
% subplot(2,2,2)
% plot(grating_range,TOD(1,:))
% subplot(2,2,3)
% plot(prism_range,GDD(2,:))
% subplot(2,2,4)
% plot(prism_range,TOD(2,:))

% figure
% [A,H1,H2] = plotxx(grating_range,GDD(1,:),prism_range,GDD(2,:));
% ylim(A(1),[-8e4,-1e4])
% ylim(A(2),[-8e4,-1e4])
% title('GDD plot for grating and prism')
% %legend('grating','prism')
% figure
% [A,H1,H2] = plotxx(grating_range,TOD(1,:),prism_range,TOD(2,:));
% ylim(A(1),[-16e4,14e4])
% ylim(A(2),[-16e4,14e4])
% title('TOD plot for grating and prism')
% %legend('grating','prism')