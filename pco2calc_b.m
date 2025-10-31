function [pco2,co3,hco3,co2aq,pH]=pco2calc_b(DIC,Alk,T,S)
% [pco2, co3,hco3,co2aq,pH]=pco2calc_b(DIC,Alk,T,S)
% Calculates the total dissolved inorganic carbon, carbonate ion
% bicarbonate ion, aqueous CO2 and pH given atmospheric pCO2 in atm,
% alkalinity in mol/kg, temperature in C and salinity in PSU. 
% *Accounting for the impact of borate*

%Start by defining T,S



e=exp(1);

tk=T+273.15;

one=ones(size(T));



%Dissociation constant of water



lnkw=148.9802-13847.26./tk-23.6521*log(tk)+S.^(1/2).*(-5.977+118.67./tk+1.0495*log(tk))-0.01615*S;


kw=e.^(lnkw);



%Henry's Law constant

lnkh=-60.2409+93.4517*100.0./tk+23.3585*log(tk/100)+S.*(0.023517-0.023656*tk/100+0.0047036*(tk/100).^2);
kh=e.^(lnkh);




%Carbonic acid to bicarbonate



nlk1=-62.008+3670.7./tk+9.7944*log(tk)-0.0118*S+0.000116*S.^2;


k1=10.^(-nlk1);




% Bicarbonate to carbonate



nlk2=4.777+1394.7./tk-0.0184*S+0.000118*S.^2;


k2=10.^(-nlk2);



%Borate

lnkb=(1.0./tk).*(-8966.9-2890.53*S.^(1/2)-77.942*S+1.728*S.^(1.5)-0.0996*S.^2)+148.0248+137.1942*S.^(0.5)+1.62142*S+0.053105*S.^(0.5).*tk+log(tk).*(-24.4344-25.085*S.^(0.5)-0.2474*S);

kb=e.^(lnkb);

TB=11.85e-6*S;

a(:,1)=-ones(size(T));
a(:,2)=-(Alk+kb+k1);
a(:,3)=(-Alk.*(kb+k1)-kb.*k1-k1.*k2+kw+k1.*DIC+kb.*TB);
a(:,4)=(-Alk.*(k1.*k2+k1.*kb)-k1.*k2.*kb+kw.*(kb+k1)+2*k1.*k2.*DIC+k1.*kb.*DIC+k1.*kb.*TB);
a(:,5)=(-Alk.*k1.*k2.*kb+kw.*k1.*kb+kw.*k1.*k2+2*k1.*k2.*kb.*DIC+k1.*k2.*kb.*TB);
a(:,6)=kw.*kb.*k1.*k2;

for j=1:max(size(T))
  r=roots(a(j,:));
  h=max(r(r>0));
  pH(j)=-log10(h);
  hco3(j)=k1(j)*h*DIC(j)/(h^2+k1(j)*h+k1(j)*k2(j));
  co3(j)=k2(j)*hco3(j)/h;
  co2aq(j)=h*hco3(j)/k1(j);
  pco2(j)=co2aq(j)/kh(j);
end
return