function [v] = Qmultiplication(Q, b)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Funkcja pomocnicza, która mnoży macierz Q (produkt kolejnych 
% macierzy Givensa) z rozkładu QR przez wektor b
% WEJŚCIE:
%   Q - macierz wyjściowa funkcji QRfactor (zawierająca pary wartości 
%       {c_k, s_k} obrotów Givensa w kolejnych kolumnach
%   b - wektor
% WYJŚCIE:
%   v - wektor rozwiązania

n = length(b);
v = b;

for k = 1: n-1
    b(k) = Q(1,k) * v(k) - Q(2,k) * v(k+1);
    b(k+1) = Q(2,k) * v(k) + (Q(1,k)) * v(k+1);
    v = b;
end

end % function