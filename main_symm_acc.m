clear
filename='out_symm_7.txt';
fid=fopen(filename,'w');

M = 5000; N = 500;
% M = 5; N = 1;
fprintf(fid,'%s %d %s %d\n','M=',M,'N=',N);
fprintf(fid,'\n');
fprintf(fid,'%3s %3s %2s %4s %8s','set','n','p','a','power');
fprintf(fid,'\n');

C = [0,0];
n1 = [40,60,80];
p1 = [2,4,6];
setting1 = [7];
C(2) = length(n1)*length(p1)*length(setting1);

for setting = setting1
  for n = n1
    for p = p1

      % test 1
      power = test_symm1_acc(M,N,setting,n,p);
      for i = 1:4
        fprintf(fid,'%3d %3d %2d %4.1f %8.4f\n',setting,n,p,i/2,power(i));
      end

      % test 2
      power = test_symm2_acc(M,N,setting,n,p);
      for i = 1:3
        fprintf(fid,'%3d %3d %2d %4.1f %8.4f\n',setting,n,p,i/2,power(i));
      end

      C(1) = C(1) + 1

    end
  end
end

fclose(fid);
