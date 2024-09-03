# BSfate is a bistable switch circuits prediction method for directional single-cell RNA-seq data. It is written in R language (R >= 4.2.1).

**BSfate** takes a scRNA-seq dataset from the directional differentiation process as input. First, pseudo-time information of cells can be obtained through trajectory inference methods such as slingshot. BSFate models the gene expression as a nonlinear function of pseudo-time and utilizes a nonlinear least squares method to screen TFs that display switch-like or transient activation patterns. By combinatorically paired, BSFate forms candidate gene pairs and calculates the significance score for each candidate. The pairs with the top rank are anticipated to constitute the core cell fate circuits. Meanwhile, TFs can be ordered based on the most prominent rank of their respective gene pairs for further investigation or experimental validation.

!(“C:/Users/温柔/Desktop/jpg/1.jpg”)
