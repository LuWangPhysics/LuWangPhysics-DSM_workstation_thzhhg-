folder_name='DSM_THzHHG_LuWang/my_output/E_01_5.00E5V_per_m77k_1507fs_449em3_1.000gamma_dz500e15/'
z={"10"};
for num=1:length(z)
thz_spec=importdata(folder_name+z{num}+"nmthz_e_w.txt"); 
f_thz=importdata(folder_name+"thz_omega.txt")./(2*pi*1e12);
e_thz_w=thz_spec(:,1)+1i.*thz_spec(:,2);
N_half=fix(length(e_thz_w)/2);
e_thz_w(1:N_half)=flip(conj(e_thz_w((N_half+1):end)));


figure(1)
hold on 
plot(f_thz,abs(e_thz_w))
set(gca,'Yscale','Log')
xlim([0,20])
title('Output total spectra')
xlabel('f (THz)')
ylabel('|E(f)|')


% 
 figure(2)
eff=importdata(folder_name+z{num}+"nmeff_3hg.txt");
z0=importdata(folder_name+"z.txt");
 plot(z0(1:length(eff)).*1e6,eff)
 xlabel('z (um)')
 ylabel('efficiency of 3HG (%)')

 figure(3)
 e_t=importdata(folder_name+z{num}+"nmet_thz_total.txt");
 t=importdata(folder_name+"t.txt");
 plot(t.*1e12,fftshift(e_t))
 xlim([-5,5])
 xlabel('t (ps)')
 ylabel('E(t)')
 title('totoal electric field')

  %e field of the 3HG
 df=f_thz(2)-f_thz(1);
 f_3hg_range=fix(3/df)+(-fix(0.85/df):1:fix(1/df));
 ew_3hg=zeros(size(e_thz_w));
 ew_3hg(f_3hg_range+N_half)=e_thz_w(f_3hg_range+N_half);
 ew_3hg(1:N_half)=flip(conj(ew_3hg((N_half+1):end)));
 
 et_3hg=fftshift(ifft(ifftshift(ew_3hg))).*(2*N_half*df*1e12);
 figure(4)
plot(t.*1e12,real(et_3hg))
 xlim([-5,5])
 xlabel('t (ps)')
 ylabel('E(t)')
title('E(t) of the third harmonic only')
end

