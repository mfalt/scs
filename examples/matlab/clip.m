function y = clip(x,u,l)
	y = max(min(x,u),l);
end