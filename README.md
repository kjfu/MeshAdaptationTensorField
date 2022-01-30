<!--
 * @Author: Kejie Fu
 * @Date: 2022-01-08 22:29:56
 * @LastEditTime: 2022-01-08 22:48:09
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /MeshMetric/README.md
-->
# MeshMetric
Generate a metric tensor (3x3 matrix) filed from a scalar field defined on a input tetrahedral mesh, which can be used in mesh adaptation.

# Usage
## 1.Basic Usage
Input *TempName.vtk* and *TempName.sol*, **MeshMetrix** will output *TempName.o.sol*
```shell
$ MeshMetric TempName
```
## 2.Options
### Output visualize file in *.vtk
Input *TempName.vtk* and *TempName.sol*, **MeshMetrix** will output *TempName.o.sol* and *TempName.o.vtk* 
```shell
$ MeshMetrix -v TempName
```
