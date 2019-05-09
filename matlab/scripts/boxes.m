function boxes = boxes(h, w)
    

    img1 = zeros(h, w);
    img2 = zeros(h, w);
    img3 = zeros(h, w);
    img4 = zeros(h, w);
    
    h2 = h/2 + 1;
    w2 = w/2 + 1;
    h = h;
    w = w;

    img1(1:w2, 1:h2) = 1;
    img2(w2+1:w, 1:h2) = 1;
    img3(1:w2, h2+1:h) = 1;
    img4(w2+1:w, h2+1:h) = 1;
    
%     figure()
%     hold on
%     axis equal
%     imagesc(img1)
%     
%     figure()
%     hold on
%     axis equal
%     imagesc(img2)
%     
%     figure()
%     hold on
%     axis equal
%     imagesc(img3)
%     
%     figure()
%     hold on
%     axis equal
%     imagesc(img4)
    
    boxes = {img1, img2, img3, img4};
    
end

