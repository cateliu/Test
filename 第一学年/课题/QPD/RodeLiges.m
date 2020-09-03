function T=RodeLiges(theta,v1,k,o)
v = v1-o;
v_para = sum(v.*k)*k/norm(k);
T = (1-cos(theta))*v_para+cos(theta)*v+sin(theta)*cross(k,v);
end