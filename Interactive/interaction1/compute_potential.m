function Usave = compute_potential(Q,g,m,r)
%COMPUTE_POTENTIAL
%    USAVE = COMPUTE_POTENTIAL(Q,G,M,R)

%    This function was generated by the Symbolic Math Toolbox version 7.2.
%    24-Jan-2018 22:52:17

Usave = g.*m.*r.*sin(Q);
