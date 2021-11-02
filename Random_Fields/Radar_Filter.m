function [Rr] = Radar_Filter(r, phi, x_, y_, x, y, dx, dy, Rt, Nt)


for ti = 1:Nt
    for i = 1:length(r)
        for l = 1:length(phi)
            F = squeeze(Rt(ti, :, :));
            if (i == length(r)) && (l == length(phi))
                xil1(i, l) = r(i) .* cos(phi(l));
                yil1(i, l) = r(i) .* sin(phi(l));
                xil2(i, l) = r(i) .* cos(phi(l));
                yil2(i, l) = r(i) .* sin(phi(l));

                [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
                [~, idy1(i, l)] = min(abs(y_ - yil1(i, l)));
                [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
                [~, idy2(i, l)] = min(abs(y_ - yil2(i, l)));

                idxmin(i, l) = min(idx1(i, l), idx2(i, l));
                idxmax(i, l) = max(idx1(i, l), idx2(i, l));

                idymin(i, l) = min(idy1(i, l), idy2(i, l));
                idymax(i, l) = max(idy1(i, l), idy2(i, l));

                Rr(ti, i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;

                A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));

                D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                    (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);

            elseif (i == length(r)) && (l ~= length(phi))
                xil1(i, l) = r(i) .* cos(phi(l));
                yil1(i, l) = r(i) .* sin(phi(l));
                xil2(i, l) = r(i) .* cos(phi(l + 1));
                yil2(i, l) = r(i) .* sin(phi(l + 1));

                [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
                [~, idy1(i, l)] = min(abs(y_ - yil1(i, l)));
                [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
                [~, idy2(i, l)] = min(abs(y_ - yil2(i, l)));

                idxmin(i, l) = min(idx1(i, l), idx2(i, l));
                idxmax(i, l) = max(idx1(i, l), idx2(i, l));

                idymin(i, l) = min(idy1(i, l), idy2(i, l));
                idymax(i, l) = max(idy1(i, l), idy2(i, l));

                Rr(ti, i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;

                A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));

                D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                    (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);


            elseif (i ~= length(r)) && (l == length(phi))
                xil1(i, l) = r(i) .* cos(phi(l));
                yil1(i, l) = r(i) .* sin(phi(l));
                xil2(i, l) = r(i + 1) .* cos(phi(l));
                yil2(i, l) = r(i + 1) .* sin(phi(l));

                [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
                [~, idy1(i, l)] = min(abs(y_ - yil1(i, l)));
                [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
                [~, idy2(i, l)] = min(abs(y_ - yil2(i, l)));

                idxmin(i, l) = min(idx1(i, l), idx2(i, l));
                idxmax(i, l) = max(idx1(i, l), idx2(i, l));

                idymin(i, l) = min(idy1(i, l), idy2(i, l));
                idymax(i, l) = max(idy1(i, l), idy2(i, l));

                Rr(ti, i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;
               
                A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));

                D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                    (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);

            else

                xil1(i, l) = r(i) .* cos(phi(l));
                yil1(i, l) = r(i) .* sin(phi(l));
                xil2(i, l) = r(i + 1) .* cos(phi(l + 1));
                yil2(i, l) = r(i + 1) .* sin(phi(l + 1));

                [~, idx1(i, l)] = min(abs(x_ - xil1(i, l)));
                [~, idy1(i, l)] = min(abs(y_ - yil1(i, l)));
                [~, idx2(i, l)] = min(abs(x_ - xil2(i, l)));
                [~, idy2(i, l)] = min(abs(y_ - yil2(i, l)));

                idxmin(i, l) = min(idx1(i, l), idx2(i, l));
                idxmax(i, l) = max(idx1(i, l), idx2(i, l));

                idymin(i, l) = min(idy1(i, l), idy2(i, l));
                idymax(i, l) = max(idy1(i, l), idy2(i, l));

                Rr(ti, i, l) = sum(sum(abs(F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l))))) .* dx .* dy;

                A(i, l).a(:, :) = F(idxmin(i, l):idxmax(i, l), idymin(i, l):idymax(i, l));

                D(i, l).d(:, :) = sqrt((x(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2 + ...
                    (y(idymin(i, l):idymax(i, l), idxmin(i, l):idxmax(i, l))).^2);

            end
        end
    end
    Rr(ti, :, :) = squeeze(Rr(ti, :, :))./max(max(Rr(ti, :, :)));
end
