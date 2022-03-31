222

111
%вариант Шефера 28_05_21
%бинарные сигналы
%глазковая диаграмма
%прямоугольный импульс
%эллиптический фильтр 
%вычисляется во временной области импульс, выходящий из фильтра, 
%затем вычисляется сигнал
%можно учесть дисперсию: вычисляется преобразование Фурье сигнала,
%оно умножается на передаточную функцию канала, затем вычисляется 
%обратное преобразование Фурье
clear all;
close all;
set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
%-------------------------------------------------------------------
N1=16;
N=16384*64;%число точек преобразования Фурье
NN=N/N1;%число битов
L_trasa=40; %50; 100; в километрах

EbN0 = 16;
D=17;% в пс/км/нм - дисперсионный волноводный параметр
lambda=1550;% в нм
v=3*10^8;%электродинамическая постоянная
alpha=pi*D*lambda^2*L_trasa*10^(-21)/v;%
x=randi([0,1],NN,1);

switch 1    % переключение формата кодирования

    case 1
         b = zeros(N,1);
         d = x;
         d_d = not(d);

         for i = 1:NN
             b(i + 1) = xor (d_d(i),b(i)); % прекодер
             c(i) = b(i) + b(i + 1) - 1;   % двубинарный код
         end
% % ------------------------------------------
         data_1 = repmat(c',1,N1)';
    case 2
         data_1 = repmat(x,1,N1)';
    case 3
         ww(1)=1;   
         for i = 1:NN
             ww(i+1)=ww(i)+x(i);
             if ww(i+1)==2
                ww(i+1)=0;
             end
         end   
         for i = 1:NN
             b(i)=ww(i)+ww(i+1);
             c(i)=b(i)-1;
         end
         data_1 = repmat(c',1,N1)';
%          data_1 = repmat(ww',1,N1)';
end

% % ------------пропускание прямоугольного импульса через фильтр Бесселя
T=1e-10;

wmax=(pi/T)*N1;

dt=T/N1;


data_2 = reshape(data_1,size(data_1,1)*size(data_1,2),1); %первичный цифровой сигнал
data_1 = repmat(x,1,N1)';
data_2_bin = reshape(data_1,size(data_1,1)*size(data_1,2),1); %первичный цифровой сигнал
mm = 0;
jj = 1;
% s_filtr(N:1) = zeros;



% for long=L_trasa:-1:1
%     long
%     L = long;
%     mm = mm + 1;
% %     jj = jj + 1;
%     long_koef = 0.07:-(0.04 / (L_trasa)):0.03;
% 
%     [b,a] = ellip(1,0.5,60,0.07);
%     [b_1,a_1] = ellip(1,0.5,5,long_koef(mm));
% 
%     filtr_on = 1;
% 
%     if filtr_on
%         data_2 = filtfilt(b,a,data_2);
%         data_3 = filtfilt(b_1,a_1,data_2);
%         data_2_bin = filtfilt(b,a,data_2_bin);
%         data_3_bin = filtfilt(b_1,a_1,data_2_bin);
%     else
%         data_3 = data_2;
%         data_3_bin = data_2_bin;  
%     end
%     data_2_orig = data_2;
% 
%     zatuhanie = 1.02*L;
% 
%     % for snr_ber = -10:1:40
%     % mm=mm+1;
%     % EbN0 = snr_ber;
%     snr=EbN0-10*log10(0.5*T/dt);
%   
%     data_3_shum = awgn(data_3,snr,'measured');
%     data_3_bin_shum = awgn(data_3_bin,snr,'measured');
%     data_2 = data_3_shum;
%     data_3_shum = data_3_shum / zatuhanie;
%     data_3_bin_shum = data_3_bin_shum / zatuhanie;
% 
%     % k_CD = 10; % Коэф уширения
%     % m_CD_1 = 1:1/k_CD:((N1/4)*1/k_CD - 1/k_CD) + 1;
%     % m_CD_2 = 1:-1/k_CD:-((N1/4)*1/k_CD - 1/k_CD) + 1;
%     % m_CD_3 = fliplr(m_CD_2);
%     % m_CD_4 = fliplr(m_CD_1);
%     % mm_CD = [m_CD_1, m_CD_2, m_CD_3, m_CD_4];
%     % mm_CD = mm_CD';
%     % mm_CD
%     % data_test = zeros(N,0);
%     % data_test(1:16) = mm_CD;
%     % for i = 1: (N / N1) - 1
%     %     data_test(i * N1+1:i * N1 + 16) = mm_CD; 
%     % end
%     % data_test = data_test';
%     % data_test_3 = data_2 .* data_test;
%     % data_test = repmat(mm_CD',1,2)';
%     % data_test_2 = reshape(data_test,size(data_test,1)*size(data_test,2),1); %первичный цифровой сигнал
% 
% 
%     S=fft(data_2_orig,N);
%     S1(1:length(S)/2)=S(length(S)/2+1:length(S));
%     S1(length(S)/2+1:length(S))=S(1:length(S)/2); 
%     M=max(abs(S1));
%     P_orig=20*log10(abs(S1)/M);
% 
%     S=fft(data_3_shum,N);
%     S1(1:length(S)/2)=S(length(S)/2+1:length(S));
%     S1(length(S)/2+1:length(S))=S(1:length(S)/2); 
%     M=max(abs(S1));
%     P_chanel=20*log10(abs(S1)/M);
% 
% 
%     dw=2*wmax/N;
%     w=0:dw:wmax-dw;%+
%     w1=-wmax:dw:-dw;%-
%     w2=[w1, w];
%     w2 = w2/ 10000000000;
% 
% 
%     SS=S';
% 
%     S1(1:length(S)/2)=SS(length(S)/2+1:length(S));
%     S1(length(S)/2+1:length(S))=SS(1:length(S)/2); 
%     M=max(abs(S1));
%     P=20*log10(abs(S1)/M);
% 
%     s=ifft(SS',N);
%     s_filtr = s;
%     s_bin_filtr = data_3_bin_shum;
%     if filtr_on
%         s_filtr = filtfilt(b,a,s);
%         s_bin_filtr = filtfilt(b,a,data_3_bin_shum);
%     end 
% 
% 
%     %%%%%%%%%%%%%%%%%----------ПРИЕММММ----------%%%%%%%%%%%%%%
%     Gamma = 0.3*max(s_filtr);
%     Gamma_bin = 0.3*max(s_bin_filtr);
%     count = 0;
%     j = 1;
%     s_dmd(1:NN) = zeros;
%     s_dmd_bin(1:NN) = zeros;
% 
%     for i = 1:N
%         count = count + 1;
%         s_dmd(j) = s_dmd(j) + s_filtr(i);
%         s_dmd_bin(j) = s_dmd_bin(j) + s_bin_filtr(i);
%         if count == N1
%             count = 0;
%             s_dmd(j) = s_dmd(j)/N1;
%             s_dmd_bin(j) = s_dmd_bin(j)/N1;
%             if s_dmd(j) > Gamma
%                 s_dmd(j)= 1; 
%             elseif s_dmd(j) < -Gamma
%                 s_dmd(j)= -1;     
%             else
%                 s_dmd(j)= 0;
%             end
%             if s_dmd_bin(j) > Gamma_bin
%                 s_dmd_bin(j)= 1;    
%             else
%                 s_dmd_bin(j)= 0;
%             end
%             j = j + 1;
%         end 
%     end
%     % s_dmd = s_dmd ./N1; 
%     c = abs (c) ;
%     s_dmd = abs (s_dmd);
%     [er,ratio] = biterr(c, s_dmd); % эта функция сравнивает b0 
%     [er_bin,ratio_bin] = biterr(c, s_dmd_bin); % эта функция сравнивает b0 
%     %c b и выводит в  er число несовпадающих бит
% 
%     ber1(mm,jj)=er/NN;%-10*log10(er/N);
%     ber1_bin(mm,jj)=er_bin/NN;%-10*log10(er/N);  
%     Snr(mm,jj)= L_trasa - long;
% end

[b,a] = ellip(1,0.5,60,0.07);
[b_1,a_1] = ellip(1,0.5,5,0.07);

filtr_on = 1;

if filtr_on
   data_2 = filtfilt(b,a,data_2);
   data_3 = filtfilt(b_1,a_1,data_2);
   data_2_bin = filtfilt(b,a,data_2_bin);
   data_3_bin = filtfilt(b_1,a_1,data_2_bin);
else
   data_3 = data_2;
   data_3_bin = data_2_bin;  
end
data_2_orig = data_2;

zatuhanie = 1.02*2;

snr=EbN0-10*log10(0.5*T/dt);
  
data_3_shum = awgn(data_3,snr,'measured');
data_3_bin_shum = awgn(data_3_bin,snr,'measured');
data_2 = data_3_shum;
data_3_shum = data_3_shum / zatuhanie;
data_3_bin_shum = data_3_bin_shum / zatuhanie;


s_filtr = data_3_shum;
s_bin_filtr = data_3_bin_shum;
if filtr_on
   s_filtr = filtfilt(b,a,data_3_shum);
   s_bin_filtr = filtfilt(b,a,data_3_bin_shum);
end 



















% % Snr = Snr .* -1;
% ber1        = ber1 / 10;
% ber1_bin    = ber1_bin / 10;
% semilogy(Snr,ber1,'kx')
% hold on;
% semilogy(Snr,ber1_bin,'ko')
% legend({'Бинарный сигнал','Двубинарный сигнал'}, 'location', 'southeast')
% xlabel('Протяженность ЛС, км')
% ylabel('BER')



kol_vo_bit_plot = 7;        %Число бит для отрисовки
t_plot = 1:1/N1:kol_vo_bit_plot + 1-1/N1;
f_plot = -N:1/N1:N;

if 1    % осциллограф
    figure('position',[20 680 750 300])
    hold on
    plot(t_plot, data_2_orig(1:kol_vo_bit_plot * N1),'k--','linewidth',2);
    plot(t_plot, s_filtr(1:kol_vo_bit_plot * N1),'k.','linewidth',2);
    plot(t_plot, data_3_shum(1:kol_vo_bit_plot * N1),'k','linewidth',2);
    grid
    xlabel('Информ. последовательность,биты')
    ylabel('Амплитуда')
    legend({'Первичный цифровой сигнал','Приянятый и отфтльтрованный сигнал','Сигнал пройденный через канал связи'}, 'location', 'southeast')
end

% fmax=1/2/dt;
% wmax=2*pi*fmax;
% dw=wmax/N;
% w0=1/T/4;
% w1=0:dw:wmax-dw;

if 0    % Спектр
    figure('position',[1130 680 750 300]);
    hold on;
    subplot(1,2,1);
    plot(w2,P_orig,'k','linewidth',0.1);
    grid;
    axis([-10 10 -80 0]);
    xlabel('Частота, ГГц');
    ylabel('Амплитуда, дБ');
    subplot(1,2,2);
    plot(w2,P_chanel,'k','linewidth',0.1);
    grid;
    axis([-10 10 -80 0]);
    xlabel('Частота, ГГц');
    ylabel('Амплитуда, дБ');
end
if 0    % Глазковая Диаграмма
%     s_filtr = s_filtr - 1;
    H = eyediagram(abs(s_filtr),2*N1,2*T,0,'k');
    xlabel('Время, сек.');
    ylabel('Амплитуда');
end

