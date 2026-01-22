def filter_sparse_cloud(chunk, rec_uncert=10, proj_acc=3,
                        reproj_err=0.5, remove_percent=10,
                        min_cluster_size=30):
    f = Metashape.PointCloud.Filter()
    f.init(chunk, Metashape.PointCloud.Filter.ReconstructionUncertainty)
    f.removePoints(threshold=rec_uncert)
    f.init(chunk, Metashape.PointCloud.Filter.ProjectionAccuracy)
    f.removePoints(threshold=proj_acc)
    f.init(chunk, Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(threshold=reproj_err)
    f.init(chunk, Metashape.PointCloud.Filter.ReprojectionError)
    f.removePoints(percentage=remove_percent)
    chunk.point_cloud.removeSmallComponents(size=min_cluster_size)

