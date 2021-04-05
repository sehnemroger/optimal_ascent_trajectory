clc; clear all;
j=1;
for i=0:.1:2*pi
u=2; 
beta=i
v=3; we=4; inc=3.5; r=10;

alfa1=atan2((-u*cos(beta)+(v-r*we*cos(inc))*sin(beta)),((v-r*we*cos(inc))*cos(beta)+u*sin(beta)));
alfa2=beta-atan2(u,(v-r*we*cos(inc)));
result=alfa1-alfa2
end