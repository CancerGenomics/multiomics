compare.rows.cnv <- function(calculated.wcgna, calculated.legacy, row) {
  expect_equal(calculated.wcgna[row,1], calculated.legacy[row,1])
  # Disabled as getGeneLocation has beed modified to return "unk"
  #expect_equal(calculated.wcgna[row,2], calculated.legacy[row,2])
  expect_equal(round(as.numeric(calculated.wcgna[row,3]), digits=3), round(as.numeric(calculated.legacy[row,3]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,4]), digits=3), round(as.numeric(calculated.legacy[row,4]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,5]), digits=3), round(as.numeric(calculated.legacy[row,5]), digits = 3))
}

compare.rows.mirna <- function(calculated.wcgna, calculated.legacy, row) {
  expect_equal(calculated.wcgna[row,1], calculated.legacy[row,1])
  expect_equal(calculated.wcgna[row,2], calculated.legacy[row,2])
  expect_equal(round(as.numeric(calculated.wcgna[row,3]), digits=3), round(as.numeric(calculated.legacy[row,3]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,4]), digits=3), round(as.numeric(calculated.legacy[row,4]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,5]), digits=3), round(as.numeric(calculated.legacy[row,5]), digits = 3))
}

compare.rows.meth <- function(calculated.wcgna, calculated.legacy, row) {
  expect_equal(calculated.wcgna[row,1], calculated.legacy[row,1])
  # Disabled as getGeneLocation has beed modified to return "unk"
  #expect_equal(calculated.wcgna[row,2], calculated.legacy[row,2])
  expect_equal(calculated.wcgna[row,3], calculated.legacy[row,3])
  expect_equal(round(as.numeric(calculated.wcgna[row,4]), digits=3), round(as.numeric(calculated.legacy[row,4]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,5]), digits=3), round(as.numeric(calculated.legacy[row,5]), digits = 3))
  expect_equal(round(as.numeric(calculated.wcgna[row,6]), digits=3), round(as.numeric(calculated.legacy[row,6]), digits = 3))
}