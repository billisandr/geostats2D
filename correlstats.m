function a = correlstats(x,y)
a(1,1) = mean(abs(y - x));
a(1,2) = max(abs(y - x));
a(1,3) = mean((y - x).^2);
a(1,4) = sqrt(mean((y - x).^2));
a(1,5) = corr(x,y);
a(1,6) = corr(x,y,'type','Spearman');
end

    