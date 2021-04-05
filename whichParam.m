function param=whichParam
%
% Função que escolhe o parametro a ser utilizado com entradas do usuário.
% --------------------------------------
%| Puramente organizacional |
% --------------------------------------

disp('Selecione o valor correspondente ao parametro desejado :');
op=input('(1) -> r0\n(2) -> sigma0\n(3) -> v0\n(4) -> T1\n(5) -> T2\n(6) -> T3\n(7) -> m12\n(8) -> m23\n(9) -> s12\n(10) -> s23\n(11) -> k1\n(12) -> k2\n(13) -> k3\n(14) -> rf\n(15) -> m0\n(16) -> rf\n');
switch op
    case 1
        param= 'r0';
    case 2
        param= 'sigma0';
    case 3
        param= 'v0';
    case 4
        param='T1';
    case 5
        param='T2';
    case 6
        param='T3';
    case 7
        param='m12';
    case 8
        param='m23';
    case 9
        param='s12';
    case 10
        param='s23';
    case 11
        param='k1';
    case 12
        param='k2';
    case 13
        param='k3';
    case 14
        param='rf';
    case 15
        param='m0';
    case 16
        param='vf';
end