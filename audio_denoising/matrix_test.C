void matrix_test(){
  Int_t n = 3; // n x n matrix
  Double_t a[] = {1,2,3};
  TVectorD aa;
  aa.Use(n,a); 

  TMatrixD A(n,n);
  TMatrixDRow(A,0) = aa;

  //here is a for to reallocate the elements and fill the matrix 
  for(Int_t i = 1; i<n; i++){ // this for to fill the matrix
    Double_t temp = aa[0];
    aa[0] = aa[n-1];
    for(Int_t j = 1; j<n;j++){
      Double_t temp2 = aa[j];
      aa[j] = temp; 
      temp = temp2; 
    }
    TMatrixDRow(A,i) = aa;
  }
  for(Int_t i = 0; i<n; i++){
    for(Int_t j = 0; j<n; j++){
      cout << A[i][j] << " ";
    }
    cout << "\n";
  }

}
