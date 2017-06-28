function x = dirac(N,i,a)
%DIRAC generiert einen Dirac-Impuls
% Aufruf: x = dirac(N,i,a)
% x = Dirac-Impuls
% N = Anzahl werte
% i = Position des Impulses in x
% a = Amplitude des Impulses
% Beispiel:
%             x = dirac(64,32,1)
%
x(1,1:N) = zeros(1,N) ;
x(i) = a ;
