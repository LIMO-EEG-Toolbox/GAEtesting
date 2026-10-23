function [L, num] = bwlabeln(BW, conn)
% 2D 4-connected labels. Stand-in when IPT / SPM mex are missing.
% Enough for LIMO Time clustering (limo_findcluster labels each channel's time axis).
if nargin < 2
    conn = 4;
end
if ndims(BW) > 2 && ~ismatrix(BW)
    error('bwlabeln shim is 2D only (LIMO Time maps). Need IPT or SPM mex for TF.')
end
BW = squeeze(BW) ~= 0;
[nr, nc] = size(BW);
L = zeros(nr, nc);
num = 0;
dirs = [-1 0; 1 0; 0 -1; 0 1];
for i = 1:nr
    for j = 1:nc
        if ~BW(i, j) || L(i, j) ~= 0
            continue
        end
        num = num + 1;
        qr = i;
        qc = j;
        L(i, j) = num;
        head = 1;
        while head <= numel(qr)
            r = qr(head);
            c = qc(head);
            head = head + 1;
            for d = 1:4
                rr = r + dirs(d, 1);
                cc = c + dirs(d, 2);
                if rr >= 1 && rr <= nr && cc >= 1 && cc <= nc && BW(rr, cc) && L(rr, cc) == 0
                    L(rr, cc) = num;
                    qr(end + 1) = rr; %#ok<AGROW>
                    qc(end + 1) = cc; %#ok<AGROW>
                end
            end
        end
    end
end
end
