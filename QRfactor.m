function [Q, R] = QRfactor(a, b, c)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Funkcja pomocnicza, która rozkłada macierz trójdiagonalną A na iloczyn 
% macierzy Q i R
% WEJŚCIE:
%   a - wektor elementów podprzekątnej macierzy A
%   b - wektor elementów głównej przekątnej macierzy A
%   c - wektor elementów nadprzekątnej macierzy A
% WYJŚCIE:
%   Q - macierz, której kolumnami są kolejne pary wartości {c_k, s_k}
%       obrotów Givensa
%   R - macierz, której wierszami są 3 wektory, pierwszy wiersz to wektor
%       elementów leżących na głównej przekątnej macierzy R w 
%       rozkładzie QR, drugi wiersz to elementy z nadprzekątnej, trzeci
%       wiersz to elementy z przekątnej nad nadprzekątną.

% Inicjalizacja zmiennych pomocniczych i parametrów wyjściowych
n = length(b);
d = zeros(1, n-2);
R = zeros(3, n);
Q = zeros(2, n-1);

for k = 1: n-1
    % Obliczenie macierzy Givensa
    [c_k, s_k] = givens(b(k), a(k));
    G = [c_k -s_k; s_k c_k];
    Q(1, k) = c_k;
    Q(2, k) = s_k;
    
    % Zastosowanie macierzy Givensa do obrotu macierzy A
    m = min([k+2, n]);
    if m == k+2
        T = [b(k), c(k), 0; a(k), b(k+1), c(k+1)];
        T = G * T;
        c(k+1) = T(2, 3);
        d(k) = T(1, 3);
    else
        T = [b(k), c(k); a(k), b(k+1)];
        T = G * T;
    end
    a(k) = T(2, 1);
    b(k) = T(1, 1);
    c(k) = T(1, 2);
    b(k+1) = T(2, 2);

end

% Utworzenie wektorów diagonali macierzy R
R(1, :) = b;
R(2, 1:n-1) = c;
R(3, 1:n-2) = d;

end % function