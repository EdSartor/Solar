close all

%----------------------------------------
%---------------USER INPUTS--------------
%----------------------------------------
filename = 'astmg173.csv';
P = 5; %Weighting factor for Integrated Power Match
I = 1; %Weighting factor for Irradiance Match
Q = 1; %0 to integrate from infrared, 1 to integrate from UV

%Recommended Parameters P = 5, I = 1, Q = 1



% ----------------Import LED Spectra---------
A = importdata('LEDspectra2.dat');
lam = A(:,1);
LED = zeros(149,24);

%------------------Plot LED Spectra-----------
figure(1)
hold on
for i = 1:24
    L(:,i) = A(:,i+1);
    plot(lam,L(:,i));
end
title('LED Spectra')
xlabel('Wavelength (nm)')
ylabel('Spectral Irradiance (W m-2 nm-1)')
hold off

%-----------------Import and Condense Desired Spectrum------
data = importdata(filename);
dlam = 5;

data = importdata(filename);
len = length(data(:,1));
B = zeros(149,2);
B(:,1) = lam;
j = 1;
for i = 1:149
    while (j <= len & data(j,1) <= lam(i))
        if round(data(j,1),0) == A(i,1)
            B(i,2) = data(j,2);
            j = j + 1;
        else j = j + 1;
        end
    end
end

%Add extra power from below 360 nm (LED cutoff point) to the spectrum
k = find(data(:,1)==360);
UVpow = trapz(data(1:k,2))*(data(2,1)-data(1,1));
B(1,2) = B(1,2)+UVpow;

%----------------Generate Integrated Spectrum----------
intB = zeros(149,2);
%intB = zeros(149,2);
intB(:,1) = lam;
for i = 1:149
    intB(i,2) = trapz(B(1:i,2))*dlam;
end

%---------------Generate Reverse Integrated Spectrum-
intBr = zeros(149,2);
for i = 1:149
    intBr(i,2) = trapz(B(149-i+1:149,2))*dlam;
    intBr(i,1) = lam(149-i+1);
end



%-------Fit Spectra----------------
x0 = ones(1,24);
UB = [70,67,47,68,66,100,96,42,83,141,103,34,105,24,115,94,42,50,165,1,2,0,0,0];
LB = ones(1,24)*0;
Param = [P,I,Q];
x2 = fmincon(@(x)fitfunc(x,A,B,Param),x0,[],[],[],[],LB,UB);

x2 = round(x2,0);
%-------Plot Spectra, Fit, and Difference-------
C = x2*L';
C = C';

diff = abs(C-B(:,2));

intC = zeros(149,2);
intC(:,1) = lam;
for i = 1:149
    intC(i,2) = trapz(C(1:i,1))*dlam;
end

intCr = zeros(149,2); 
for i = 1:149
    intCr(i,2) = trapz(C(149-i+1:149,1))*dlam;
    intCr(i,1) = lam(149-i+1);
end

figure(3)
hold on
plot(lam(2:149),B(2:149,2))
plot(lam,C(:,1))
plot(lam(2:149),diff(2:149,1))
title('Spectrum and Fit')
xlabel('Wavelength (nm)')
ylabel('Spectral Irradiance (W m-2 nm-1)')
legend('Spectrum','Fit','Difference')
hold off

figure(4)
hold on
if Q == 1
plot(lam,intC(:,2))
plot(lam,intB(:,2))
else
end
if Q == 0
plot(intBr(:,1),intCr(:,2))
plot(intBr(:,1),intBr(:,2))
else
end
title('Integrated Power Comparison')
xlabel('Wavelength (nm)')
ylabel('Power (W m-2)')
legend('Fit','Spectrum')
hold off

n = 1:24;
A = [n; x2]

fileID = fopen('SpectraMax.dat','w');
fprintf(fileID,'%6s %12s\n','Spectrum:',filename);
fprintf(fileID,'%6s %12s\n','LED #','Value');
fprintf(fileID,'%6.2f %12.8f\n',A);
fclose(fileID);

function cost = fitfunc(x,A,B,Param)
Q = Param(3);
lam = A(:,1);
leng=length(lam);
L = zeros(leng,24);
for i = 1:24
    L(:,i) = A(:,i+1);
end

intB = zeros(leng,2);  %---------INTEGRATE SPECTRUM
intB(:,1) = lam;
for i = 1:leng
    intB(i,2) = trapz(B(1:i,2));
end

intBr = zeros(leng,2);  %-------------REVERSE INTEGRATE SPECTRUM
for i = 1:leng
    intBr(i,2) = trapz(B(leng-i+1:leng,2));
    intBr(i,1) = lam(leng-i+1);
end

C = x*L'; %--------FIT LEDS
C = C';

intC = zeros(leng,2); %----------INTEGRATE FIT
intC(:,1) = lam;
for i = 1:leng
    intC(i,2) = trapz(C(1:i,1));
end

intCr = zeros(leng,2);  %-------------REVERSE INTEGRATE FIT
for i = 1:leng
    intCr(i,2) = trapz(C(leng-i+1:leng,1));
    intCr(i,1) = lam(leng-i+1);
end

diff = abs(C-B(:,2));   %------------------CALCULATE COST FUNCTION
intdiff = abs(intB(:,2)-intC(:,2));
revintdiff = abs(intBr(:,2)-intCr(:,2));
cost = Param(2)*sum(diff)+Param(1)*(Q*sum(intdiff)-(Q-1)*sum(revintdiff));
end