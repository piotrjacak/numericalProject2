function [c, s] = givens(a, b)
% Projekt 2, zadanie 13
% Piotr Jacak, 327354

% Funkcja pomocnicza do obliczania pary [c, s] wartości obrotów Givensa
% WEJŚCIE:
%   a, b - dwie wartości rzeczywiste (pewne dwa elementy macierzy)
% WYJŚCIE:
%   c, s - dwie wartości (c = cos(teta), s = sin(teta)) obrotów Givensa

if b == 0
    c = 1;
    s = 0;
else
    if abs(b) > abs(a)
        t = -a / b;
        s = 1 / sqrt(1 + t^2);
        c = s * t;
    else
        t = -b / a;
        c = 1 / sqrt(1 + t^2);
        s = c * t;
    end
end

end % function