function imout = imresize_3d(imin,newsize)

if size(imin,2)~=newsize(2)
    imout1 = zeros(newsize(1),newsize(2),size(imin,3));
    for n=1:size(imin,3)
        imout1(:,:,n) = imresize(imin(:,:,n),[newsize(1) newsize(2)]);
    end
else
    imout1=imin;
end

imout = zeros(newsize(1),newsize(3),newsize(2));

imout1 = permute(imout1,[1 3 2]);

for n = 1:newsize(2)

    imout(:,:,n) = imresize(imout1(:,:,n),[newsize(1) newsize(3)]);

end

imout = permute(imout,[1 3 2]);


end