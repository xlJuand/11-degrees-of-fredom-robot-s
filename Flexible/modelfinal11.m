%% 11 degrees of fredom robot's
clear
clc
% Model variables
% Relative joint angles
syms th1 th2 th3 th4 th5 th6 th7 th8 th9 tht thp
q=[th1;th2;th3;th4;th5;th6;th7;th8;th9;tht;thp];
%Relative joint velocities
syms dth1 dth2 dth3 dth4 dth5 dth6 dth7 dth8 dth9 dtht dthp
dq=[dth1;dth2;dth3;dth4;dth5;dth6;dth7;dth8;dth9;dtht;dthp];
%Gravity
syms g
% link lengths
syms L1 L2 L3 L4 L5 L6 L7 L8 L9 Lt Lp
% link masses
syms m1 m2 m3 m4 m5 m6 m7 m8 m9 mt mp
%link inertias
syms Izz_1 Izz_2 Izz_3 Izz_4 Izz_5 Izz_6 Izz_7 Izz_8 Izz_9 Izz_t Izz_p
%% Kinetic energy
% K[i] = 0.5*m[i]*v[ci]^2 + 0.5*I[i]*dq_i^2
% V[ci] : linear velocity
%important points
p1 = 2*L1*[cos(th1);sin(th1)];
p2 = p1 + 2*L2*[cos(th2);sin(th2)];
p3 = p2 + 2*L3*[cos(th3);sin(th3)];
p4 = p3 + 2*L4*[cos(th4);sin(th4)];
p5 = p4 + 2*L5*[cos(th5);sin(th5)];
p6 = p5 + 2*L6*[cos(th6);sin(th6)];
p7 = p6 + 2*L7*[cos(th7);sin(th7)];
p8 = p7 + 2*L8*[cos(th8);sin(th8)];
p9 = p8 + 2*L9*[cos(th9);sin(th9)];
pt = p9 + 2*Lt*[cos(tht);sin(tht)];
pp = pt + 2*Lp*[cos(thp);sin(thp)];

%Position of center of masses
r1 = L1*[cos(th1);sin(th1)];
r2 = p1 + L2*[cos(th2);sin(th2)];
r3 = p2 + L3*[cos(th3);sin(th3)];
r4 = p3 + L4*[cos(th4);sin(th4)];
r5 = p4 + L5*[cos(th5);sin(th5)];
r6 = p5 + L6*[cos(th6);sin(th6)];
r7 = p6 + L7*[cos(th7);sin(th7)];
r8 = p7 + L8*[cos(th8);sin(th8)];
r9 = p8 + L9*[cos(th9);sin(th9)];
rt = p9 + Lt*[cos(tht);sin(tht)];
rp = pt + Lp*[cos(thp);sin(thp)];


%Center of masses velocities
dr1=jacobian(r1,q)*dq;
dr2=jacobian(r2,q)*dq;
dr3=jacobian(r3,q)*dq;
dr4=jacobian(r4,q)*dq;
dr5=jacobian(r5,q)*dq;
dr6=jacobian(r6,q)*dq;
dr7=jacobian(r7,q)*dq;
dr8=jacobian(r8,q)*dq;
dr9=jacobian(r9,q)*dq;
drt=jacobian(rt,q)*dq;
drp=jacobian(rp,q)*dq;

%% Kinetic energy of links
KE_m1=0.5*m1*(dr1.'*dr1)+0.5*Izz_1*(dth1)^2;
KE_m1=simplify(KE_m1);

KE_m2=0.5*m2*(dr2.'*dr2)+0.5*Izz_2*(dth2)^2;
KE_m2=simplify(KE_m2);

KE_m3=0.5*m3*(dr3.'*dr3)+0.5*Izz_3*(dth3)^2;
KE_m3=simplify(KE_m3);

KE_m4=0.5*m4*(dr4.'*dr4)+0.5*Izz_4*(dth4)^2;
KE_m4=simplify(KE_m4);

KE_m5=0.5*m5*(dr5.'*dr5)+0.5*Izz_5*(dth5)^2;
KE_m5=simplify(KE_m5);


KE_m6=0.5*m6*(dr6.'*dr6)+0.5*Izz_6*(dth6)^2;
KE_m6=simplify(KE_m6);

KE_m7=0.5*m7*(dr7.'*dr7)+0.5*Izz_7*(dth7)^2;
KE_m7=simplify(KE_m7);

KE_m8=0.5*m8*(dr8.'*dr8)+0.5*Izz_8*(dth8)^2;
KE_m8=simplify(KE_m8);

KE_m9=0.5*m9*(dr9.'*dr9)+0.5*Izz_9*(dth9)^2;
KE_m9=simplify(KE_m9);

KE_mt=0.5*mt*(drt.'*drt)+0.5*Izz_t*(dtht)^2;
KE_mt=simplify(KE_mt);

KE_mp=0.5*mp*(drt.'*drt)+0.5*Izz_p*(dthp)^2;
KE_mp=simplify(KE_mp);

%total kinetic energy
KE=KE_m1+KE_m2+KE_m3+KE_m4+KE_m5+KE_m6+KE_m7+KE_m8+KE_m9+KE_mt+KE_mp;
KE=simplify(KE);
%% Potential energy
PE1 = m1*g*(L1+r1(2));
PE1=simplify(PE1);

PE2 =m2*g*(2*L1+L2+r2(2));
PE2=simplify(PE2);

PE3 =m3*g*(2*L1+2*L2+L3+r3(2));
PE3=simplify(PE3);

PE4 =m4*g*(2*L1+2*L2+2*L3+L4+r4(2));
PE4=simplify(PE4);

PE5 =m5*g*(2*L1+2*L2+2*L3+2*L4+L5+r5(2));
PE5=simplify(PE5);

PE6 =m6*g*(2*L1+2*L2+2*L3+2*L4+2*L5+L6+r6(2));
PE6=simplify(PE6);

PE7 =m7*g*(2*L1+2*L2+2*L3+2*L4+2*L5+2*L6+L7+r7(2));
PE7=simplify(PE7);

PE8 =m8*g*(2*L1+2*L2+2*L3+2*L4+2*L5+2*L6+2*L7+L8+r8(2));
PE8=simplify(PE8);

PE9 =m9*g*(2*L1+2*L2+2*L3+2*L4+2*L5+2*L6+2*L7+2*L8+L9+r9(2));
PE9=simplify(PE9);

PEt=mt*g*(2*L1+2*L2+2*L3+2*L4+2*L5+2*L6+2*L7+2*L8+2*L9+Lt+rt(2));
PEt=simplify(PEt);

PEp=mp*g*(2*L1+2*L2+2*L3+2*L4+2*L5+2*L6+2*L7+2*L8+2*L9+2*Lt+2*rt+Lp+rp(2));
PEp=simplify(PEp);

PE=PE1+PE2+PE3+PE4+PE5+PE6+PE7+PE8+PE9+PEt+PEp;
PE=simplify(PE);

%model matrices
G=jacobian(PE,q).';
G=simplify(G);

%mass-inertial matrix
D_mtx = simplify(jacobian(KE,dq).');
D_mtx = simplify(jacobian(D_mtx,dq));
%Coriolis and centrifugal matrix
syms C_mtx real
n=max(size(q));

for k=1:n
for j=1:n
C_mtx(k,j)=0*g;
for i=1:n
C_mtx(k,j)=C_mtx(k,j)+1/2*(diff(D_mtx(k,j),q(i))+...
diff(D_mtx(k,i),q(j)) -...
diff(D_mtx(i,j),q(k)))*dq(i);
end
end
end

C_mtx=simplify(C_mtx);
%% Create Matlab functions for each EOM
param_list = { 'g','p(1)'; 'L1','p(2)'; 'L2','p(3)'; 'L3','p(4)';
'L4','p(5)'; 'L5','p(6)'; 'L6','p(7)'; 'L7','p(8)';
'L8','p(9)'; 'L9','p(10)';'Lt','p(11)';'Lp','p(12)';
'm1','p(13)';'m2','p(14)';'m3','p(15)';'m4','p(16)';
'm5','p(17)';'m6','p(18)';'m7','p(19)';'m8','p(20)';
'm9','p(21)';'mt','p(22)';'mp','p(23)';'Izz_1','p(24)';'Izz_2','p(25)';
'Izz_3','p(26)';'Izz_4','p(27)';'Izz_5','p(28)';
'Izz_6','p(29)';'Izz_7','p(30)';'Izz_8','p(31)';
'Izz_9','p(32)';'Izz_t','p(33)';'Izz_p','p(34)';
};

list_q = {'th1','q(1)'; 'th2','q(2)';'th3','q(3)'; 'th4','q(4)';'th5','q(5)'; 'th6','q(6)';'th7','q(7)'; 'th8','q(8)';'th9','q(9)';'tht','q(t)';'thp','q(p)'};
list_dq = {'dth1','q(1)'; 'dth2','q(2)';'dth3','q(3)'; 'dth4','q(4)';'dth5','q(5)'; 'dth6','q(6)';'dth7','q(7)'; 'dth8','q(8)';'dth9','q(9)';'dtht','q(t1)';'dthp','q(t2)'};
%% Model matrices

write_fcn('D_matrix.m',{'q','p'},[list_q;param_list],{D_mtx,'D'});
write_fcn('C_matrix.m',{'q','dq','p'},[list_q;list_dq;param_list],{C_mtx,'C'});
write_fcn('G_matrix.m',{'q','p'},[list_q;param_list],{G,'G'});
write_fcn('p1_vector.m',{'q','p'},[list_q;param_list],{p1,'p1'});
write_fcn('p2_vector.m',{'q','p'},[list_q;param_list],{p2,'p2'});
write_fcn('p3_vector.m',{'q','p'},[list_q;param_list],{p3,'p3'});
write_fcn('p4_vector.m',{'q','p'},[list_q;param_list],{p4,'p4'});
write_fcn('p5_vector.m',{'q','p'},[list_q;param_list],{p5,'p5'});
write_fcn('p6_vector.m',{'q','p'},[list_q;param_list],{p6,'p6'});
write_fcn('p7_vector.m',{'q','p'},[list_q;param_list],{p7,'p7'});
write_fcn('p8_vector.m',{'q','p'},[list_q;param_list],{p8,'p8'});
write_fcn('p9_vector.m',{'q','p'},[list_q;param_list],{p9,'p9'});
write_fcn('pt_vector.m',{'q','p'},[list_q;param_list],{pt,'pt'});
write_fcn('pp_vector.m',{'q','p'},[list_q;param_list],{pp,'pp'});