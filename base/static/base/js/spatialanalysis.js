document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_spatialanalysis_results = await fetchAPI(
    "/api/preloadspatialanalysis",
    0,
    csrftoken
  );
  document.getElementById("preloadresult").textContent =
    preload_spatialanalysis_results.preload_result;
  console.log(preload_spatialanalysis_results);
});
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const spatialanalysis_results = await fetchAPI(
    "/api/spatialanalysis",
    0,
    csrftoken
  );
  console.log(spatialanalysis_results);
  toggleLoading(false, "processbutton");
  document.getElementById("imgbox").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");
}
