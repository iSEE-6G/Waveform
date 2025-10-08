P = [1 1 0; 0 0 1; 1 0 1; 0 1 0];
G = [eye(4) P];
H = [P' eye(3)];

Gnew = G;
Gnew(1,:) = xor(G(1,:), G(2, :));
Gnew(3,:) = xor(G(2,:), G(3, :));

Hnew = H;
Hnew(1,:) = xor(H(1,:), H(2, :));
Hnew(3,:) = xor(H(2,:), H(3, :));

flag = 1;
while flag
    GG = randi([0,1], 3, 7);
    tmp = mod(Hnew*GG',2);
    if tmp == eye(3)
        flag = 0;
    end
end