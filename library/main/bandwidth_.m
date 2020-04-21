function nband = bandwidth_(g)
% !
% ! This function finds the element bandwidth from g.
% !
nband = max(g(g>0)) - min(g(g>0));
end