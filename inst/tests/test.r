context('Tmem68 mediation scan LOD scores')

test_that("equals to saved.version",{

  data("Tmem68")
  med <- mediation.scan(target=Tmem68$target,
                        mediator=Tmem68$mediator,
                        annotation=Tmem68$annotation,
                        covar=Tmem68$covar,
                        qtl.geno=Tmem68$qtl.geno,
                        method="double-lod-diff")

  data("Tmem68.lod")

  expect_equal(med$LOD, Tmem68.lod)
})
