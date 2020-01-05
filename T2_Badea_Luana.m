clc %clc curata fereastra de mesaje la fiecare noua rulare a programului
T=40;
f=1/T;
w=2*pi*f;
t = 0:0.0001:2*T;
%voi reprezenta doua perioade din semnal, pasul de discretizare al
%semnalului continuu este dt=10^-4.
N=25;
C = zeros(1,2*N+1);
%am initializat vectorul de coeficienti cu valori nule
x=sawtooth(w*t, 0.5)
%semnal triunghiular simetric monoredresat de perioada pricipala T=32 s si
%durata D=8 s, semnal obtinut dintr-un semnal triunghiular simetric ce a
%fost monoredresat
for n = -N:N
    C(n+N+1) = 1/T * integral(@(t)(sawtooth(w*t,0.5)).*exp(-1j*n*w*t),0,T) ;
    %calculul propriu-zis al coeficientiilor cu formula analitica
    re = real(C(n+N+1));
    im = imag(C(n+N+1));
    if abs(re)<10^-10
        re = 0;
     end
    if abs(im)<10^-10
        im = 0;
    end
    %daca partea reala sau imaginara a unui coeficient este extrem de mica
    %atunci o voi neglija
    C(n+N+1) = re+1j*im ;
    %ca index intr-un vector nu pot avea valori negative de aceea indexul
    %fiecarui element din vector este cu 51 mai mare ca indexul teoretic ai
    %coeficientului, exemplu:C(-50)|teoretic=C(1)|in program, insa acest
    %lucru nu modifica valorile coeficientiilor si functionalitatea
    %programului
end
xr = 0;
for n = -N:N
    xr = xr + C(n+N+1)*exp(1j*n*w*t) ;
end
%am recostruit semnalul initial cu 100 de componente
figure(1);
hold on
plot(t,xr);%reprezentarea semnalului reconstruit
plot(t,x);%reprezentarea semnalului initial
xlabel('Timpul[s]');
ylabel('(t)-linie continua si xr(t)-linie intrerupra');
title('Suprapunerea semnalelor x(t) si xr(t)')
%semnalul triunghiular se poate aproxima cu mai putine componente( in jur
%de 15-20) si astfel diferenta dintre x si xr este aproape inexistenta
hold off
%prin hold on/off pot reprezenta mai multe grafice in acelasi sistem de
%coordonate
figure(2);
hold on
stem((-N:N)*w,2*abs(C));
%functia stem este utilizata pentru reprezentarea functilor sau a
%seturi de date cu valori discrete
plot((-N:N)*w,2*abs(C),'-go');%infasuratoarea realizate din segmente de dreapta
xlabel('Pulsatia w [rad/s]');
ylabel('Amplitudinile Ak=2|C(nw)|');

title('Spectrul de Amplitudini');
hold off
%In concluzie, analiza Fourier a semnalelor analogice(continue) ne permite
%sa exprimam orice semnal ce indeplineste criteriul de suficienta Diriclet
%intr-o suma finita de semnale elementare lucru ce este folositor in
%analiza comoda a circuitelor in domeniul fazorial sau reconstruirea unui
%semnal necunoscut pe baza spectrului sau de amplitudini si de faze ce pot
%fi aflate cu un analizator de spectru