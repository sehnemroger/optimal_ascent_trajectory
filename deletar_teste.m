%deletar_teste
clc; clear all; close all;
erro=1e-10:1e-11:2e-5;
er_lim=1e-8; % erro ref;
frac_er=(er_lim-erro)/er_lim;
y3 =-1+1*2*sigmf(frac_er,[.003 0]);
%plot(erro_corrigido,linspace(0,length(erro),length(erro)),'LineWidth',3)
plot(frac_er,y3,'LineWidth',3)