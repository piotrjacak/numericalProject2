function [eigenValue, errest, i] = ...
    P2Z13_PJA_findEigenvalue(a, b, c, u, tol, maxIter)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Program znajduje przybliżenie wartości własnej trójdiagonalnej 
% macierzy rzeczywistej A, leżącej najbliżej podanej wartości u 
% metodą odwróconej metody potęgowej z normowaniem. 
% Do rozwiązywania układów równań używa rozkładu macierzy A-u*I 
% wyznaczonego za pomocą rotacji Givensa. 
% Program kończy działanie po osiągnięciu zadanej dokładności lub
% w przypadku przekroczenia maksymalnej liczby iteracji
% WEJŚCIE:
%   a - wektor elementów rzeczywistych leżących na podprzekątnej 
%       macierzy A (długość wektora należy do przedziału [0, 199999]) 
%   b - wektor elementów rzeczywistych leżących na głównej przekątnej 
%       macierzy A (długość wektora należy do przedziału [1, 200000]) 
%   c - wektor elementów rzeczywistych leżących na nadprzekątnej 
%       macierzy A (długość wektora należy do przedziału [0, 199999])
%   u - wartość rzeczywista leżąca najbliżej szukanej wartości własnej
%   tol - tolerancja na błąd występująca w warunku zakończenia obliczeń
%   maxIter - maksymalna liczba iteracji metody (na wypadek braku
%             zbieżności metody potęgowej)
% WYJŚCIE:
%   eigenValue - znaleziona wartość własna
%   errest - oszacowanie błędu w warunku zakończenia obliczeń
%   i - końcowa liczba wykonanych iteracji

n = length(b);
% Inicjalizacja wektora początkowego x elementami losowymi
x = randn(n, 1);
errest = Inf;
i = 1;

% Rozkład QR macierzy A - u*I 
[Q, R] = QRfactor(a, b - u, c);

while i < maxIter && errest > tol
    % Normalizacja wektora x
    xnorm = x./norm(x);
    % Obliczenie kolejnego przybliżenia wektora własnego
    % korzystając z odwrotnej metody potęgowej
    v = Qmultiplication(Q, xnorm);
    x = upperTriangularSolve(R(1, :), R(2, :), R(3, :), v);
    % Obliczenie wartości własnej na podstawie wektora własnego
    eigenValue = (conj(x)' * multiplication(a, b-u, c, x)) ...
        / (conj(x)' * x) + u;
    % Obliczenie błędu przybliżenia
    errest = norm(x - eigenValue .* xnorm) / abs(eigenValue);
    i = i + 1;
end

end % function