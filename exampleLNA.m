clear all
close all


%%% X ->{k1} 2X
%%% X + Y ->{k2} 2Y
%%% Y ->{k3} 

%%% X ->{k1} c1*lambda1 (where c1 =0,1,2)
%%% lambda1 + lambda2 ->{k1} c2*Y
%%% lambda2 ->{k3}  


syms k1 k2 k3 X Y covX covY covXY real

alpha1=k1*X
alpha2=k2*X*Y
alpha3=k3*Y


%%%Species are encoded as [X Y]

v1=[1 0]

v2=[-1 1]

v3=[0 -1]



F=v1*alpha1+v2*alpha2+v3*alpha3




J=jacobian(F,[X Y]) 

G=v1'*v1*alpha1
G=G+(v2'*v2*alpha2)
G=G+(v3'*v3*alpha3)

C=[covX covXY ;
   covXY covY ]


dCovdt=J*C+C*(J')+G %Differential equation COvariance

dCovdt=simplify(dCovdt)