function f = fuison_base_parts(b1, b2)
% average strategy

w1 = 0.5;
w2 = 0.5;
f = w1.*b1+w2.*b2;
% figure;imshow(f);
end
