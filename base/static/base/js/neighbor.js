document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_neighbor_result = await fetchAPI(
    "/api/preloadneighbor",
    0,
    csrftoken
  );
  console.log(preload_neighbor_result);
  document.getElementById("preloadresult").textContent =
    preload_neighbor_result.preload_result;
});
