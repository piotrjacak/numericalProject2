function [x] = multiplication(a, b, c, v)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Funkcja pomocnicza, która mnoży macierz trójdiagonalną A
% przez wektor v
% WEJŚCIE:
%   a - wektor elementów podprzekątnej macierzy A
%   b - wektor elementów głównej przekątnej macierzy A
%   c - wektor elementów nadprzekątnej macierzy A
%   v - wektor o długości równej długości wektora b
% WYJŚCIE:
%   x - wektor rozwiązania

n = length(b);
x = zeros(n, 1);

x(1) = b(1) * v(1) + c(1) * v(2);
for k = 2: n-1
    x(k) = a(k-1) * v(k-1) + b(k) * v(k) + c(k) * v(k+1);
end
x(n) = a(n-1) * v(n-1) + b(n) * v(n);

end % function