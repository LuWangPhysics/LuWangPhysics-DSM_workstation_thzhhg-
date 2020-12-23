z={"0","270"};
for num=1:length(z)
ir_1=importdata(z{num}+"nmir1_ef_abs.txt");
ir_2=importdata(z{num}+"nmir2_ef_abs.txt");
thz_spec=importdata(z{num}+"nmthz_ef_abs.txt");
energy_cons=importdata("energy_cons.txt"); 
eff=importdata(z{num}+"nmeff.txt");

t=importdata("t.txt");
f_ir1=importdata("ir1_omega.txt")./(2*pi*1e12);
f_ir2=importdata("ir2_omega.txt")./(2*pi*1e12);
f_thz=importdata("thz_omega.txt")./(2*pi*1e12);

figure(1)
hold on
plot(f_ir1,sqrt(ir_1(:,1).^2+ir_1(:,2).^2))
xlabel('f (THz)')
ylabel('$|E(f)|$')
figure(2)
hold on
plot(f_ir2,sqrt(ir_2(:,1).^2+ir_2(:,2).^2))
xlabel('f (THz)')
ylabel('$|E(f)|$')
figure(3)
hold on 
plot(f_thz,sqrt(thz_spec(:,1).^2+thz_spec(:,2).^2))
xlim([0,20])
xlabel('f (THz)')
ylabel('$|E(f)|$')


% 
% f_thz=omega_thz./(2*pi*1e12);
% J_thz=importdata("J_thz0nm.txt");
% J_thz1=importdata("J_thz10nm.txt");
% J_ir1=importdata("J_ir10nm.txt");
% J_ir2=importdata("J_ir20nm.txt");
% J_ir11=importdata("J_ir110nm.txt");
% J_ir21=importdata("J_ir210nm.txt");
% rm_1=0;
% J_ir1_3=J_ir1-rm_1*J_ir11;
% J_ir2_3=J_ir2-rm_1*J_ir21;
% J_thz_3=J_thz-rm_1*J_thz1;
% 
% subplot(1,3,1); 
% hold on
% plot(f_ir1,sqrt((J_ir1_3(:,1).^2+J_ir1_3(:,2).^2)))
% set(gca,'Yscale','Log')
% xlabel('f (THz)')
% ylabel('|J(\omega)| (C/m^2)')
% subplot(1,3,2); 
% hold on
% plot(f_ir2,sqrt((J_ir2_3(:,1).^2+J_ir2_3(:,2).^2)))
% set(gca,'Yscale','Log')
% xlabel('f (THz)')
% ylabel('|J(\omega)| (C/m^2)')
% subplot(1,3,3); 
% hold on
% plot(f_thz,sqrt((J_thz(:,1).^2+J_thz(:,2).^2)))
% xlabel('f (THz)')
% ylabel('|J(\omega)| (C/m^2)')
% df=(f_ir1(2)-f_ir1(1))*1e12;
% J_ir1t=fftshift(ifft(ifftshift(J_ir1(:,1)+J_ir1(:,2)))).*length(J_ir1)*df;
% J_ir2t=fftshift(ifft(ifftshift(J_ir2(:,1)+J_ir2(:,2)))).*length(J_ir2)*df;
% figure
% plot(t.*1e15,J_ir1t+J_ir2t)
% 
% figure(4)
% eff=importdata(z{num}+"nmeff.txt");
% z0=importdata("z.txt");
% plot(z0.*1e6,eff)
% xlabel('z (um)')
% 
% figure(5)
% e_t=importdata(z{num}+"nmet_thz.txt");
% t=importdata("t.txt");
% plot(t.*1e12,fftshift(e_t))
% xlim([-5,5])
% xlabel('t (ps)')
% ylabel('E(t)')

end

figure(5)
hold on
plot(eff)
figure(6)
hold on
plot(energy_cons)