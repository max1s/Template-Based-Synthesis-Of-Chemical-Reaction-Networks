clear all
close all


%%% X ->{k1} 2X
%%% X + Y ->{k2} 2Y
%%% Y ->{k3} 

%%% X ->{k1} c1*lambda1 (where c1 =0,1,2)
%%% lambda1 + lambda2 ->{k1} c2*Y
%%% lambda2 ->{k3}  


syms k1 k2 k3 X Y covX covY real

alpha1=k1*X
alpha2=k2*X*Y
alpha3=k3*Y


%%%Species are encoded as [X Y]

v1=[1 0]

v2=[-1 1]

v3=[0 -1]



F=v1*alpha1+v2*alpha2+v3*alpha3

pause;


%F1=(c1+c2*lambda1A)*alpha1 +    (-1+c4*lambda3A)*alpha2
%F2=c2*lambda1B*alpha1 + (c4*lambda3B)*alpha2


%F=[-alpha1+alpha2;-alpha2+alpha1-vX2*alpha1] % drift matrix
J=jacobian(F,[X]) 
pause;
G=v1'*v1*alpha1
G=G+(v2'*v2*alpha2)
G=G+(v3'*v3*alpha3)

C=[covX covY ;
   covY covX ]

pause;


dCovdt=J*C+C*(J')+G %Differential equation COvariance

dCovdt=simplify(dCovdt)



