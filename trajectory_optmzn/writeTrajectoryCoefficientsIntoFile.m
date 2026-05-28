function writeTrajectoryCoefficientsIntoFile(a, b, c_pol)
% -----------------------------------------------------------------------
% Write trajectory coefficients (Fourier + 5th-order polynomial) to scripts.
% Supports arbitrary number of joints (7-DOF for XM7p).
% -----------------------------------------------------------------------
n = size(a, 1);
prcsn = 10;

a_coefs = cell(n, 1);
b_coefs = cell(n, 1);
c_coefs = cell(n, 1);
for j = 1:n
    a_coefs{j} = sprintf('a%d = [', j);
    b_coefs{j} = sprintf('b%d = [', j);
    c_coefs{j} = sprintf('c%d = [', j);
end

for i = 1:size(a, 2)
  if i < size(a, 2)
    suffix = ',';
  else
    suffix = ']\n';
  end
  for j = 1:n
    a_coefs{j} = strcat(a_coefs{j}, num2str(a(j, i), prcsn), suffix);
    b_coefs{j} = strcat(b_coefs{j}, num2str(b(j, i), prcsn), suffix);
  end
end

for i = 1:size(c_pol, 2)
  if i < size(c_pol, 2)
    suffix = ',';
  else
    suffix = ']\n';
  end
  for j = 1:n
    c_coefs{j} = strcat(c_coefs{j}, num2str(c_pol(j, i), prcsn + 2), suffix);
  end
end

fileID_a = fopen('trajectory_optmzn/coeffs4_UR/a_coeffs.script', 'w');
fileID_b = fopen('trajectory_optmzn/coeffs4_UR/b_coeffs.script', 'w');
fileID_c = fopen('trajectory_optmzn/coeffs4_UR/c_coeffs.script', 'w');
for i = 1:n
    fprintf(fileID_a, a_coefs{i});
    fprintf(fileID_b, b_coefs{i});
    fprintf(fileID_c, c_coefs{i});
end
fclose(fileID_a);
fclose(fileID_b);
fclose(fileID_c);
