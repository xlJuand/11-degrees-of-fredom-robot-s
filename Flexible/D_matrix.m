function [D] = D_matrix(q,p)

  D(1,1)=p(24) + p(2)^2*p(13) + 4*p(2)^2*p(14) + 4*p(2)^2*p(15) + 4*p(2)^2*p(16) + 4*p(2)^2*p(17) +...
          4*p(2)^2*p(18) + 4*p(2)^2*p(19) + 4*p(2)^2*p(20) + 4*p(2)^2*p(21) + 4*p(2)^2*p(23) + 4*p(2)^2*p(22);
  D(1,2)=2*p(2)*p(3)*cos(q(1) - q(2))*(p(14) + 2*p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(1,3)=2*p(2)*p(4)*cos(q(1) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  D(1,4)=2*p(2)*p(5)*cos(q(1) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(1,5)=2*p(2)*p(6)*cos(q(1) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(1,6)=2*p(2)*p(7)*cos(q(1) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(1,7)=2*p(2)*p(8)*cos(q(1) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(1,8)=2*p(2)*p(9)*cos(q(1) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(1,9)=2*p(2)*p(10)*cos(q(1) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(1,10)=2*p(2)*p(11)*cos(q(1) - q(t))*(p(23) + p(22));
  D(1,11)=0;
  D(2,1)=2*p(2)*p(3)*cos(q(1) - q(2))*(p(14) + 2*p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*...
         p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(2,2)=p(25) + p(3)^2*p(14) + 4*p(3)^2*p(15) + 4*p(3)^2*p(16) + 4*p(3)^2*p(17) + 4*p(3)^2*p(18) +...
          4*p(3)^2*p(19) + 4*p(3)^2*p(20) + 4*p(3)^2*p(21) + 4*p(3)^2*p(23) + 4*p(3)^2*p(22);
  D(2,3)=2*p(3)*p(4)*cos(q(2) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  D(2,4)=2*p(3)*p(5)*cos(q(2) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(2,5)=2*p(3)*p(6)*cos(q(2) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(2,6)=2*p(3)*p(7)*cos(q(2) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(2,7)=2*p(3)*p(8)*cos(q(2) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(2,8)=2*p(3)*p(9)*cos(q(2) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(2,9)=2*p(3)*p(10)*cos(q(2) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(2,10)=2*p(3)*p(11)*cos(q(2) - q(t))*(p(23) + p(22));
  D(2,11)=0;
  D(3,1)=2*p(2)*p(4)*cos(q(1) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  D(3,2)=2*p(3)*p(4)*cos(q(2) - q(3))*(p(15) + 2*p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*...
         p(21) + 2*p(23) + 2*p(22));
  D(3,3)=p(26) + p(4)^2*p(15) + 4*p(4)^2*p(16) + 4*p(4)^2*p(17) + 4*p(4)^2*p(18) + 4*p(4)^2*p(19) +...
          4*p(4)^2*p(20) + 4*p(4)^2*p(21) + 4*p(4)^2*p(23) + 4*p(4)^2*p(22);
  D(3,4)=2*p(4)*p(5)*cos(q(3) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(3,5)=2*p(4)*p(6)*cos(q(3) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(3,6)=2*p(4)*p(7)*cos(q(3) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(3,7)=2*p(4)*p(8)*cos(q(3) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(3,8)=2*p(4)*p(9)*cos(q(3) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(3,9)=2*p(4)*p(10)*cos(q(3) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(3,10)=2*p(4)*p(11)*cos(q(3) - q(t))*(p(23) + p(22));
  D(3,11)=0;
  D(4,1)=2*p(2)*p(5)*cos(q(1) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(4,2)=2*p(3)*p(5)*cos(q(2) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(4,3)=2*p(4)*p(5)*cos(q(3) - q(4))*(p(16) + 2*p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*...
         p(23) + 2*p(22));
  D(4,4)=p(27) + p(5)^2*p(16) + 4*p(5)^2*p(17) + 4*p(5)^2*p(18) + 4*p(5)^2*p(19) + 4*p(5)^2*p(20) +...
          4*p(5)^2*p(21) + 4*p(5)^2*p(23) + 4*p(5)^2*p(22);
  D(4,5)=2*p(5)*p(6)*cos(q(4) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(4,6)=2*p(5)*p(7)*cos(q(4) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(4,7)=2*p(5)*p(8)*cos(q(4) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(4,8)=2*p(5)*p(9)*cos(q(4) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(4,9)=2*p(5)*p(10)*cos(q(4) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(4,10)=2*p(5)*p(11)*cos(q(4) - q(t))*(p(23) + p(22));
  D(4,11)=0;
  D(5,1)=2*p(2)*p(6)*cos(q(1) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,2)=2*p(3)*p(6)*cos(q(2) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,3)=2*p(4)*p(6)*cos(q(3) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,4)=2*p(5)*p(6)*cos(q(4) - q(5))*(p(17) + 2*p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,5)=p(28) + p(6)^2*p(17) + 4*p(6)^2*p(18) + 4*p(6)^2*p(19) + 4*p(6)^2*p(20) + 4*p(6)^2*p(21) +...
          4*p(6)^2*p(23) + 4*p(6)^2*p(22);
  D(5,6)=2*p(6)*p(7)*cos(q(5) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,7)=2*p(6)*p(8)*cos(q(5) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,8)=2*p(6)*p(9)*cos(q(5) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(5,9)=2*p(6)*p(10)*cos(q(5) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(5,10)=2*p(6)*p(11)*cos(q(5) - q(t))*(p(23) + p(22));
  D(5,11)=0;
  D(6,1)=2*p(2)*p(7)*cos(q(1) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,2)=2*p(3)*p(7)*cos(q(2) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,3)=2*p(4)*p(7)*cos(q(3) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,4)=2*p(5)*p(7)*cos(q(4) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,5)=2*p(6)*p(7)*cos(q(5) - q(6))*(p(18) + 2*p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,6)=p(29) + p(7)^2*p(18) + 4*p(7)^2*p(19) + 4*p(7)^2*p(20) + 4*p(7)^2*p(21) + 4*p(7)^2*p(23) +...
          4*p(7)^2*p(22);
  D(6,7)=2*p(7)*p(8)*cos(q(6) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,8)=2*p(7)*p(9)*cos(q(6) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(6,9)=2*p(7)*p(10)*cos(q(6) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(6,10)=2*p(7)*p(11)*cos(q(6) - q(t))*(p(23) + p(22));
  D(6,11)=0;
  D(7,1)=2*p(2)*p(8)*cos(q(1) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,2)=2*p(3)*p(8)*cos(q(2) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,3)=2*p(4)*p(8)*cos(q(3) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,4)=2*p(5)*p(8)*cos(q(4) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,5)=2*p(6)*p(8)*cos(q(5) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,6)=2*p(7)*p(8)*cos(q(6) - q(7))*(p(19) + 2*p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,7)=p(30) + p(8)^2*p(19) + 4*p(8)^2*p(20) + 4*p(8)^2*p(21) + 4*p(8)^2*p(23) + 4*p(8)^2*p(22);
  D(7,8)=2*p(8)*p(9)*cos(q(7) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(7,9)=2*p(8)*p(10)*cos(q(7) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(7,10)=2*p(8)*p(11)*cos(q(7) - q(t))*(p(23) + p(22));
  D(7,11)=0;
  D(8,1)=2*p(2)*p(9)*cos(q(1) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,2)=2*p(3)*p(9)*cos(q(2) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,3)=2*p(4)*p(9)*cos(q(3) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,4)=2*p(5)*p(9)*cos(q(4) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,5)=2*p(6)*p(9)*cos(q(5) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,6)=2*p(7)*p(9)*cos(q(6) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,7)=2*p(8)*p(9)*cos(q(7) - q(8))*(p(20) + 2*p(21) + 2*p(23) + 2*p(22));
  D(8,8)=p(31) + p(9)^2*p(20) + 4*p(9)^2*p(21) + 4*p(9)^2*p(23) + 4*p(9)^2*p(22);
  D(8,9)=2*p(9)*p(10)*cos(q(8) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(8,10)=2*p(9)*p(11)*cos(q(8) - q(t))*(p(23) + p(22));
  D(8,11)=0;
  D(9,1)=2*p(2)*p(10)*cos(q(1) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,2)=2*p(3)*p(10)*cos(q(2) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,3)=2*p(4)*p(10)*cos(q(3) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,4)=2*p(5)*p(10)*cos(q(4) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,5)=2*p(6)*p(10)*cos(q(5) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,6)=2*p(7)*p(10)*cos(q(6) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,7)=2*p(8)*p(10)*cos(q(7) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,8)=2*p(9)*p(10)*cos(q(8) - q(9))*(p(21) + 2*p(23) + 2*p(22));
  D(9,9)=p(32) + p(10)^2*p(21) + 4*p(10)^2*p(23) + 4*p(10)^2*p(22);
  D(9,10)=2*p(10)*p(11)*cos(q(9) - q(t))*(p(23) + p(22));
  D(9,11)=0;
  D(10,1)=2*p(2)*p(11)*cos(q(1) - q(t))*(p(23) + p(22));
  D(10,2)=2*p(3)*p(11)*cos(q(2) - q(t))*(p(23) + p(22));
  D(10,3)=2*p(4)*p(11)*cos(q(3) - q(t))*(p(23) + p(22));
  D(10,4)=2*p(5)*p(11)*cos(q(4) - q(t))*(p(23) + p(22));
  D(10,5)=2*p(6)*p(11)*cos(q(5) - q(t))*(p(23) + p(22));
  D(10,6)=2*p(7)*p(11)*cos(q(6) - q(t))*(p(23) + p(22));
  D(10,7)=2*p(8)*p(11)*cos(q(7) - q(t))*(p(23) + p(22));
  D(10,8)=2*p(9)*p(11)*cos(q(8) - q(t))*(p(23) + p(22));
  D(10,9)=2*p(10)*p(11)*cos(q(9) - q(t))*(p(23) + p(22));
  D(10,10)=p(33) + p(11)^2*p(23) + p(11)^2*p(22);
  D(10,11)=0;
  D(11,1)=0;
  D(11,2)=0;
  D(11,3)=0;
  D(11,4)=0;
  D(11,5)=0;
  D(11,6)=0;
  D(11,7)=0;
  D(11,8)=0;
  D(11,9)=0;
  D(11,10)=0;
  D(11,11)=p(34);

 