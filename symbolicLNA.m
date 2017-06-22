function symbLNA = symbolicLNA(inspecies, inreactants, inproducts, inrates)


%%% X ->{k1} 2X
%%% X + Y ->{k2} 2Y
%%% Y ->{k3} 

%%% X ->{k1} c1*lambda1 (where c1 =0,1,2)
%%% lambda1 + lambda2 ->{k1} c2*Y
%%% lambda2 ->{k3}  


species = [];
reactants = [];
products = [];
rates = [];
propensity = [];

species

for m = 1:size(1, inspecies)
	species(1, m) = sym(inspecies(1, m));
end

end
