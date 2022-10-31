function y=renyi_entro(DATA,q)

y=log(sum(DATA.^q))/(1-q);
end