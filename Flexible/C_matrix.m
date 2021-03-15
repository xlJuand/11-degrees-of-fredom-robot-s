function [C] = C_matrix(q,dq,p)

  C(1,1)=0;
  C(1,2)=2*p(2)*p(3)*dq(2)*sin(q(1) - q(2))*(p(14) + 2*p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*...
         p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(1,3)=2*p(2)*p(4)*dq(3)*sin(q(1) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(1,4)=2*p(2)*p(5)*dq(4)*sin(q(1) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(1,5)=2*p(2)*p(6)*dq(5)*sin(q(1) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(1,6)=2*p(2)*p(7)*dq(6)*sin(q(1) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(1,7)=2*p(2)*p(8)*dq(7)*sin(q(1) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(1,8)=2*p(2)*p(9)*dq(8)*sin(q(1) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(1,9)=2*p(2)*p(10)*dq(9)*sin(q(1) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(1,10)=2*p(2)*p(11)*dq(t)*sin(q(1) - q(t))*(p(23) + p(22));
  C(1,11)=0;
  C(2,1)=-2*p(2)*p(3)*dq(1)*sin(q(1) - q(2))*(p(14) + 2*p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*...
         p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(2,2)=0;
  C(2,3)=2*p(3)*p(4)*dq(3)*sin(q(2) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(2,4)=2*p(3)*p(5)*dq(4)*sin(q(2) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(2,5)=2*p(3)*p(6)*dq(5)*sin(q(2) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(2,6)=2*p(3)*p(7)*dq(6)*sin(q(2) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(2,7)=2*p(3)*p(8)*dq(7)*sin(q(2) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(2,8)=2*p(3)*p(9)*dq(8)*sin(q(2) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(2,9)=2*p(3)*p(10)*dq(9)*sin(q(2) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(2,10)=2*p(3)*p(11)*dq(t)*sin(q(2) - q(t))*(p(23) + p(22));
  C(2,11)=0;
  C(3,1)=-2*p(2)*p(4)*dq(1)*sin(q(1) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(3,2)=-2*p(3)*p(4)*dq(2)*sin(q(2) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(3,3)=0;
  C(3,4)=2*p(4)*p(5)*dq(4)*sin(q(3) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(3,5)=2*p(4)*p(6)*dq(5)*sin(q(3) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(3,6)=2*p(4)*p(7)*dq(6)*sin(q(3) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(3,7)=2*p(4)*p(8)*dq(7)*sin(q(3) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(3,8)=2*p(4)*p(9)*dq(8)*sin(q(3) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(3,9)=2*p(4)*p(10)*dq(9)*sin(q(3) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(3,10)=2*p(4)*p(11)*dq(t)*sin(q(3) - q(t))*(p(23) + p(22));
  C(3,11)=0;
  C(4,1)=-2*p(2)*p(5)*dq(1)*sin(q(1) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(4,2)=-2*p(3)*p(5)*dq(2)*sin(q(2) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(4,3)=-2*p(4)*p(5)*dq(3)*sin(q(3) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  C(4,4)=0;
  C(4,5)=2*p(5)*p(6)*dq(5)*sin(q(4) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(4,6)=2*p(5)*p(7)*dq(6)*sin(q(4) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(4,7)=2*p(5)*p(8)*dq(7)*sin(q(4) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(4,8)=2*p(5)*p(9)*dq(8)*sin(q(4) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(4,9)=2*p(5)*p(10)*dq(9)*sin(q(4) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(4,10)=2*p(5)*p(11)*dq(t)*sin(q(4) - q(t))*(p(23) + p(22));
  C(4,11)=0;
  C(5,1)=-2*p(2)*p(6)*dq(1)*sin(q(1) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(5,2)=-2*p(3)*p(6)*dq(2)*sin(q(2) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(5,3)=-2*p(4)*p(6)*dq(3)*sin(q(3) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(5,4)=-2*p(5)*p(6)*dq(4)*sin(q(4) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  C(5,5)=0;
  C(5,6)=2*p(6)*p(7)*dq(6)*sin(q(5) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(5,7)=2*p(6)*p(8)*dq(7)*sin(q(5) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(5,8)=2*p(6)*p(9)*dq(8)*sin(q(5) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(5,9)=2*p(6)*p(10)*dq(9)*sin(q(5) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(5,10)=2*p(6)*p(11)*dq(t)*sin(q(5) - q(t))*(p(23) + p(22));
  C(5,11)=0;
  C(6,1)=-2*p(2)*p(7)*dq(1)*sin(q(1) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,2)=-2*p(3)*p(7)*dq(2)*sin(q(2) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,3)=-2*p(4)*p(7)*dq(3)*sin(q(3) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,4)=-2*p(5)*p(7)*dq(4)*sin(q(4) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,5)=-2*p(6)*p(7)*dq(5)*sin(q(5) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,6)=0;
  C(6,7)=2*p(7)*p(8)*dq(7)*sin(q(6) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,8)=2*p(7)*p(9)*dq(8)*sin(q(6) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(6,9)=2*p(7)*p(10)*dq(9)*sin(q(6) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(6,10)=2*p(7)*p(11)*dq(t)*sin(q(6) - q(t))*(p(23) + p(22));
  C(6,11)=0;
  C(7,1)=-2*p(2)*p(8)*dq(1)*sin(q(1) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,2)=-2*p(3)*p(8)*dq(2)*sin(q(2) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,3)=-2*p(4)*p(8)*dq(3)*sin(q(3) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,4)=-2*p(5)*p(8)*dq(4)*sin(q(4) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,5)=-2*p(6)*p(8)*dq(5)*sin(q(5) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,6)=-2*p(7)*p(8)*dq(6)*sin(q(6) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,7)=0;
  C(7,8)=2*p(8)*p(9)*dq(8)*sin(q(7) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(7,9)=2*p(8)*p(10)*dq(9)*sin(q(7) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(7,10)=2*p(8)*p(11)*dq(t)*sin(q(7) - q(t))*(p(23) + p(22));
  C(7,11)=0;
  C(8,1)=-2*p(2)*p(9)*dq(1)*sin(q(1) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,2)=-2*p(3)*p(9)*dq(2)*sin(q(2) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,3)=-2*p(4)*p(9)*dq(3)*sin(q(3) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,4)=-2*p(5)*p(9)*dq(4)*sin(q(4) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,5)=-2*p(6)*p(9)*dq(5)*sin(q(5) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,6)=-2*p(7)*p(9)*dq(6)*sin(q(6) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,7)=-2*p(8)*p(9)*dq(7)*sin(q(7) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  C(8,8)=0;
  C(8,9)=2*p(9)*p(10)*dq(9)*sin(q(8) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(8,10)=2*p(9)*p(11)*dq(t)*sin(q(8) - q(t))*(p(23) + p(22));
  C(8,11)=0;
  C(9,1)=-2*p(2)*p(10)*dq(1)*sin(q(1) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,2)=-2*p(3)*p(10)*dq(2)*sin(q(2) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,3)=-2*p(4)*p(10)*dq(3)*sin(q(3) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,4)=-2*p(5)*p(10)*dq(4)*sin(q(4) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,5)=-2*p(6)*p(10)*dq(5)*sin(q(5) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,6)=-2*p(7)*p(10)*dq(6)*sin(q(6) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,7)=-2*p(8)*p(10)*dq(7)*sin(q(7) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,8)=-2*p(9)*p(10)*dq(8)*sin(q(8) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  C(9,9)=0;
  C(9,10)=2*p(10)*p(11)*dq(t)*sin(q(9) - q(t))*(p(23) + p(22));
  C(9,11)=0;
  C(10,1)=-2*p(2)*p(11)*dq(1)*sin(q(1) - q(t))*(p(23) + p(22));
  C(10,2)=-2*p(3)*p(11)*dq(2)*sin(q(2) - q(t))*(p(23) + p(22));
  C(10,3)=-2*p(4)*p(11)*dq(3)*sin(q(3) - q(t))*(p(23) + p(22));
  C(10,4)=-2*p(5)*p(11)*dq(4)*sin(q(4) - q(t))*(p(23) + p(22));
  C(10,5)=-2*p(6)*p(11)*dq(5)*sin(q(5) - q(t))*(p(23) + p(22));
  C(10,6)=-2*p(7)*p(11)*dq(6)*sin(q(6) - q(t))*(p(23) + p(22));
  C(10,7)=-2*p(8)*p(11)*dq(7)*sin(q(7) - q(t))*(p(23) + p(22));
  C(10,8)=-2*p(9)*p(11)*dq(8)*sin(q(8) - q(t))*(p(23) + p(22));
  C(10,9)=-2*p(10)*p(11)*dq(9)*sin(q(9) - q(t))*(p(23) + p(22));
  C(10,10)=0;
  C(10,11)=0;
  C(11,1)=0;
  C(11,2)=0;
  C(11,3)=0;
  C(11,4)=0;
  C(11,5)=0;
  C(11,6)=0;
  C(11,7)=0;
  C(11,8)=0;
  C(11,9)=0;
  C(11,10)=0;
  C(11,11)=0;

 