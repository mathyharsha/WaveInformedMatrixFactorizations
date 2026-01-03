
%max1 = max(max(WF_original(:,:,501)));
Yg = (1-sample_grid2d).*Xt;
figure(100)
for i = 1:50
   
    h=pcolor((Yg(:,:,i))); colormap('jet'); title('Original');
    h.EdgeColor = 'none'; shading interp; 
    set(gca, 'XTick', [],'YTick', []);
    %axis([1 100 1 8 min(WF_original(:)) max(WF_original(:)) -max1 max1]);
    drawnow;
   pause(0.1);
end