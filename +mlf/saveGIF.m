function saveGIF(handler,out,name)
    F         = getframe(handler); 
    [RGB,~]   = frame2im(F); 
    [IND,map] = rgb2ind(RGB,255); 
    if out == 1
        imwrite(IND,map,[name '.gif'],'gif','LoopCount',Inf); 
    else
        imwrite(IND,map,[name '.gif'],'gif','WriteMode','append','DelayTime',0); 
    end
end